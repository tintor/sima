#include <sim/system.h>

// reads Body::present, and writes to Body::future
void KinematicStep(System& system, double dtime) {
	for (unique_ptr<Body>& body : system.bodies) {
		if (body->imass == 0)
			continue;
		/*double3 acc = gGravity ? gravity : double3{0, 0, 0};
		double4 s0 = double4(body.pos.x, body.pos.y, body.vel.x, body.vel.y);
		auto s1 = RungeKutta4<dvec4, double>(s0, 0, dt, [&](dvec4 s, double t) {
			return dvec4(dvec2(s.z, s.w), acc);
		});
		body.pos = dvec2(s1.x, s1.y);
		body.vel = dvec2(s1.z, s1.w);*/
	}
}

void Apply(System& system) {
	for (auto& body : system.bodies)
		body->present = body->future;
}

inline double3 Average(cspan<double3> p) {
	double3 sum = {0, 0, 0};
	for (double3 a : p)
		sum += a;
	return sum / p.size();
}

// returns true if changes were made
bool Resolve(const Contact& contact) {
	auto& a = *contact.bodyA;
	auto& b = *contact.bodyB;

	auto& va = a.present.velocity;
	auto& vb = b.present.velocity;
	double3 n = contact.normal;

	if (dot(va, n) <= dot(vb, n))
		return false;

	auto& wa = a.present.rotation;
	auto& wb = b.present.rotation;

	double3 c = Average(contact.points); // TODO use centroid
	double3 ra = c - a.present.position;
	double3 rb = c - b.present.position;

	// velocities at contact point
	double3 vpa = va + cross(wa, ra); // TODO check order in cross product
	double3 vpb = vb + cross(wb, rb); // TODO check order in cross product

	// relative velocity at contact point
	double3 vr = vpb - vpa;

	// coefficient of restitution
	double cor = 1; // TODO comes from properties of two shapes in contact

	// reaction impulse along -normal
	double3 aim = a.imoi * cross(ra, n);
	double3 bim = b.imoi * cross(rb, n);
	double me = a.imass + b.imass + dot(cross(aim, ra) + cross(bim, rb), n);
	double jr = -(1 + cor) / me * dot(vr, n);

	// apply reaction impulse
	va -= n * (jr * a.imass);
	vb += n * (jr * b.imass);
	wa -= jr * aim;
	wb += jr * bim;

	// total external force
	double3 fe; // TODO this is needed for both bodies

	// friction tangent vector
	double3 t;
	double vrn = dot(vr, n);
	if (vrn != 0) // TODO velocity threshold needed
		t = normalize(vr - vrn * n); // TODO handle case of no tangental velocity
	else if (dot(fe, n) != 0) // TODO force threshold needed
		t = normalize(fe - dot(fe, n) * n); // TODO handle case of no tangental force
	else
		return true;

	// coefficients of friction
	double mi_s = 0.2;
	double mi_d = 0.1;

	// static / dynamic friction impulse
	double js = mi_s * jr;
	double jd = mi_d * jd;

	// friction impulse along tangent vector
	double jf;
	// TODO what is m here?
	double m;
	if (dot(vr, t) == 0 && m * dot(vr, t) <= js)
		jf = -m * dot(vr, t);
	else
		jf = -jd;

	// apply friction impulse
	va -= t * (jf * a.imass);
	vb += t * (jf * b.imass);
	wa -= jf * (a.imoi * cross(ra, t));
	wb += jf * (b.imoi * cross(rb, t));
	return true;
}

void ResolveConstraints(System& system) {
	vector<Contact> contacts;
	vector<double3> work1, work2;
	for (auto ai : range(system.bodies.size()))
		for (auto bi : range(ai + 1, system.bodies.size())) {
			Body& a = *system.bodies[ai];
			Body& b = *system.bodies[bi];

			if (a.imass + b.imass == 0)
				continue;

			// bounding circle check
			if (squared(a.present.position - b.present.position) > Tolerance * 2 + squared(a.radius + b.radius))
				continue;

			// go over each pair of convex meshes on this pair
			for (auto& sa : a.shapes)
				for (auto& sb : b.shapes) {
					Contact contact;
					double overlap;
					// TODO remember normal from last interaction of these two shapes
					int result = ClassifyConvexConvex(
						sa.mesh, sb.mesh, false, contact.normal, contact.points, overlap, work1, work2);
					if (result == 0) {
						contact.bodyA = &a;
						contact.bodyB = &b;
						contact.shapeA = &sa;
						contact.shapeB = &sb;
						contacts.push_back(std::move(contact));
					}
					if (result < 0)
						throw "objects penetrating";
				}
		}

	for (size_t iter = 0; iter < 10; iter++) {
		bool updated = false;
		for (Contact& contact : contacts)
			if (Resolve(contact))
				updated = true;
		if (!updated) {
			if (iter == 0)
				throw "resolveCollisions(): No update in the first iteration";
			break;
		}
	}
}

// finds time of first collision from present state to future state during dtime
double FindCollisionTime(System& system, double dtime) {
	double ctime = std::numeric_limits<double>::infinity();
	for (auto ai : range(system.bodies.size()))
		for (auto bi : range(ai + 1, system.bodies.size())) {
			Body& a = *system.bodies[ai];
			Body& b = *system.bodies[bi];
			// TODO
		}
	return ctime;
}

// TODO split system into groups which don't interact and advance them separately
// TODO optimize for static objects with no disturbances
void Advance(System& system, double dtime) {
	while (true) {
		// TODO compute external forces

		ResolveConstraints(system); // TODO resolve conflicting forces and velocities

		KinematicStep(system, dtime);

		double ctime = FindCollisionTime(system, dtime);
		if (ctime > dtime) {
			Apply(system);
			break;
		}

		KinematicStep(system, ctime);
		dtime -= ctime;
		Apply(system);
	}
}
