// TODO replace these with array_map and array_set
#include <unordered_map>
#include <unordered_set>

#include <geom/shape.h>
#include <core/util.h>
#include <core/timestamp.h>

#ifdef xxx
// Shape is 3d solid, immutable, purely geometric and with origin in center of mass
Shape::Shape(const Mesh3d& mesh, /*in/out*/transform3& position) {
	if (is_valid(mesh) != Validity::OK)
		throw runtime_error("invalid mesh");

	// Move to center of mass
	dvec3 center = center_of_mass(mesh);
	m_mesh = mesh; // copy
	FOR_EACH(face, m_mesh)
		face -= center;
	// TODO update 'position' out parameter, as shape is moved

	m_volume = volume(m_mesh);
	m_inertia_tensor = moment_of_inertia(m_mesh);

	// TODO HARD? rotate polyhedron, so inertia tensor only has primary axes (but return back its original orientation)
	//      so inertia tensor becomes inertia vector. Will this rotation make OBB smaller?

	m_is_convex = ::is_convex(m_mesh);

	// Bounding box
	m_is_box = ::is_box(m_mesh);
	m_box.first = min(m_mesh);
	m_box.second = max(m_mesh);

	// Compute convex edges and vertices
	std::unordered_map<segment3, dvec3> third_vertex;
	FOR_EACH(face, m_mesh)
		FOR(i, 3) {
			const dvec3& a = face[i], b = face[(i + 1) % 3], c = face[(i + 2) % 3];
			third_vertex[segment3(a, b)] = c;
		}
	std::unordered_set<dvec3> set_of_convex_vertices;
	FOR_EACH(face, m_mesh)
		FOR(i, 3) {
			const dvec3& a = face[i], b = face[(i + 1) % 3], c = face[(i + 2) % 3];
			if (edge_angle(a, b, c, third_vertex[segment3(b, a)]) < PI) {
				if (lexicographical_less(a, b))
					m_convex_edges.push_back(segment3(a, b));
			} else {
				set_of_convex_vertices.erase(a);
				set_of_convex_vertices.erase(b);
			}
		}
	FOR_EACH(v, set_of_convex_vertices)
		m_convex_vertices.push_back(v);

	// Compute radius of bounding sphere
	m_sphere_radius = 0;
	FOR_EACH(v, m_convex_vertices)
		m_sphere_radius = std::max(m_sphere_radius, squared(v));
	m_sphere_radius = sqrt(m_sphere_radius);

	// Sort faces by decreasing surface area
	std::sort(m_mesh, [](const triangle3& p, const triangle3& q) {
		return p.squared_area_x4() > q.squared_area_x4();
	});
	// Sort convex edges by decreasing length
	std::sort(m_convex_edges, [](const segment3& p, const segment3& q) {
		return squared(p.a - p.b) > squared(q.a - q.b);
	});
	// Sort convex vertices by decreasing distance from center
	std::sort(m_convex_vertices, [](const dvec3& p, const dvec3& q) {
		return squared(p) > squared(q);
	});
}

/*const auto random_directions = [](){
	std::array<dvec3, 128> dirs;
	std::default_random_engine rnd;
	FOR_EACH(d, dirs)
		d = random_direction(rnd);
	return dirs;
}();*/

real signed_distance(const dvec3& v, const Shape& shape) {
	// Can't just take the sdist sign of min_dist as two faces that share the edge can have the same dist but different sgn(sdist)
	real min_dist = std::numeric_limits<real>::max(), max_sdist = 0;
	constexpr real eps = ContactEpsilon; // ??? arbitrary, need to account for rounding errors when computing distance(vertex, face)
	FOR_EACH(face, shape.faces()) {
		plane p(face); // TODO precompute this

		real dist, sdist;
		FOR_EACH_EDGE(a, b, face)
			if (plane::sign(*a, *b, *a + p.normal, v) > 0) {
				dist = distance(v, segment3(*a, *b));
				if (dist > min_dist + eps)
					goto next;

				sdist = p.distance(v);
				goto rest;
			}

		sdist = p.distance(v);
		dist = abs(sdist);
		if (dist > min_dist + eps)
			continue;
		rest:
		if (dist < min_dist) {
			if (dist + eps < min_dist)
				max_sdist = sdist;
			min_dist = dist;
		}
		if (abs(sdist) > abs(max_sdist))
			max_sdist = sdist;
		next:;
	}
	return max_sdist < 0 ? -min_dist : min_dist;
}

real squared_distance_segment_origin(const segment3& p) {
	dvec3 d = p.b - p.a;
	real t = -dot(p.a, d);
	if (t <= 0)
		return squared(p.a);
	real dd = squared(d);
	if (t >= dd)
		return squared(p.b);
	return squared(p.a + d * (t / dd));
}

// TODO find all faces within epsilon distance from the edge
//      if any triangle is penetrated internally return true
//      if
//		rare_case:
//      collect at most one intersection point of edge with every nearby triangle
//		sort these points and check every mid point against the interior of shape
bool intersects_edge_interior(const segment3& edge, const Shape& shape) {
	// Quick bounding sphere check
	// TODO precompute right hand side and store in Shape
	if (squared_distance_segment_origin(edge) > squared(shape.sphere_radius()))
		return false;

	// Quick bounding box check
	//if (!is_edge_intersecting_box())
	//	return false;
	if (shape.is_box())
		return true;

	// Note: assumes edge vertices are outside of shape
	constexpr real eps = 1e-8; // TODO
	/*vector<triangle3> nearest;
	FOR_EACH(face, shape.faces()) {
		// if edge goes strictly through interior of triangle TODO and edge vertices are not on the plane
		if (intersects_in_point(edge, face)
				&& distance(edge, segment3(face[0], face[1])) > eps
				&& distance(edge, segment3(face[1], face[2])) > eps
				&& distance(edge, segment3(face[2], face[0])) > eps)
			return true;
		real dist = disjoint_distance(edge, face);
		if (dist < eps) {
			nearest.push_back(face);
		}
	}*/
	throw new std::runtime_error("unfinished");
}

std::pair<real, real> project_obb(dvec3 obb_position, dvec3 obb_size, dmat3 obb_orientation, dvec3 dir) {
    // TODO
	return std::pair<real, real>(0.0, 0.0);
}

bool are_oriented_boxes_intersecting(const std::pair<dvec3, dvec3>& box_a, const transform3& pos_a, const std::pair<dvec3, dvec3>& box_b, const transform3& pos_b) {
	throw new std::runtime_error("unfinished");
}

bool approx_intersects(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	// Bounding sphere check
	if (squared(pos_a.position - pos_b.position) > squared(shape_a.sphere_radius() + shape_b.sphere_radius()))
		return false;

	return are_oriented_boxes_intersecting(shape_a.box(), pos_a, shape_b.box(), pos_b);
}

bool intersects(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	// Check all vertices of one against the other
	FOR_EACH(vertex_a, shape_a.convex_vertices())
		if (intersects(pos_b.to_local(pos_a.to_global(vertex_a)), shape_b))
			return true;
	FOR_EACH(vertex_b, shape_b.convex_vertices())
		if (intersects(pos_a.to_local(pos_b.to_global(vertex_b)), shape_a))
			return true;

	// Check all edges of one against the other
	FOR_EACH(edge_a, shape_a.convex_edges())
		if (intersects_edge_interior(pos_b.to_local(pos_a.to_global(edge_a)), shape_b))
			return true;
	FOR_EACH(edge_b, shape_b.convex_edges())
		if (intersects_edge_interior(pos_a.to_local(pos_b.to_global(edge_b)), shape_a))
			return true;

	return false;
}

// TODO test intersects:
// - "simulate" two boxes "colliding" one another (different sizes and movement directions)
// - same test, but with bunnies
// - OpenGL visualization (for debugging the test)

bool is_vertex_triangle_contact(const dvec3& p, const triangle3& m, Contact& /*out*/contact) {
	dvec3 normal = Normal(m);

	// if P is outside of triangle
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + normal, p) > 0)
			return false;

	// P is inside of triangle
	dvec3 nearest = p - normal * dot(normal, p - m.a);
	contact.squared_dist = squared(nearest - p);
	if (contact.squared_dist > squared(ContactEpsilon))
		return false;
	contact.position = (nearest + p) / static_cast<real>(2);
	contact.normal = glm::normalize(normal);
	return true;
}

bool is_edge_edge_contact(const segment3& p, const segment3& q, Contact& /*out*/contact) {
	dvec3 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	real aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr real tiny = 1e-8;

	real d = aa * bb - ab * ab;
	real s = ab * bc - bb * ac;
	real t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0)) {
		dvec3 P = p.a + A * (s / d);
		dvec3 Q = q.a + B * (t / d);
		contact.squared_dist = squared(P - Q);
		if (contact.squared_dist > squared(ContactEpsilon))
			return false;
		contact.position = (P + Q) / static_cast<real>(2);
		// Note: using AxB here instead of P-Q as P and Q will be very close.
		contact.normal = glm::normalize(cross(A, B));
		if (dot(P - Q, contact.normal) < 0)
			contact.normal = -contact.normal;
		return true;
	}

	return false;
}

bool all_less_equal(const dvec3& a, const dvec3& b) {
	return a.x <= b.x && a.y <= b.y && a.z <= b.z;
}

bool approx_intersects(const dvec3& v, const Shape& shape) {
	return squared(v) <= squared(shape.sphere_radius()) && all_less_equal(shape.box().first, v) && all_less_equal(v, shape.box().second);
}

bool approx_intersects(const segment3& v, const Shape& shape) {
	// TODO
	return true;
}

// Assuming shares aren't intersecting, look for pairs of features that are within ContactEpsilon
vector<Contact> find_all_contacts(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	vector<Contact> contacts;
	Contact contact;

	// vertex vs face contacts
	FOR_EACH(vertex_a, shape_a.convex_vertices()) {
		dvec3 vertex_a_local = pos_b.to_local(pos_a.to_global(vertex_a));
		// TODO this approx check needs to be a little loose to allow for contacts
		if (!approx_intersects(vertex_a_local, shape_b))
			continue;
		FOR_EACH(face_b, shape_b.faces())
			if (is_vertex_triangle_contact(vertex_a_local, face_b, /*out*/contact)) {
				contact.position = pos_b.to_global(contact.position);
				contact.normal = pos_b.to_global_dir(contact.normal);
				contacts.push_back(contact);
			}
	}
	FOR_EACH(vertex_b, shape_b.convex_vertices()) {
		dvec3 vertex_b_local = pos_a.to_local(pos_b.to_global(vertex_b));
		// TODO this approx check needs to be a little loose to allow for contacts
		if (!approx_intersects(vertex_b_local, shape_a))
			continue;
		FOR_EACH(face_a, shape_a.faces())
			if (is_vertex_triangle_contact(vertex_b_local, face_a, /*out*/contact)) {
				contact.position = pos_a.to_global(contact.position);
				contact.normal = -pos_a.to_global_dir(contact.normal);
				contacts.push_back(contact);
			}
	}

	// edge vs edge contacts
	// TODO flip loop order if B has less edges than A
	FOR_EACH(edge_a, shape_a.convex_edges()) {
		segment3 edge_a_global = pos_a.to_global(edge_a);
		// TODO this approx check needs to be a little loose to allow for contacts
		if (!approx_intersects(pos_b.to_local(edge_a_global), shape_b))
			continue;
		// TODO avoid transforming every edge in B
		FOR_EACH(edge_b, shape_b.convex_edges())
			if (is_edge_edge_contact(edge_a_global, pos_b.to_global(edge_b), /*out*/contact))
				contacts.push_back(contact);
	}

	return contacts;
}

// returns 0 if interecting
real distance(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	// Upper bound on distance
    real dist = max(distance(pos_a.position, pos_b.position), shape_a.sphere_radius(), shape_b.sphere_radius());

	// vertex vs face
	FOR_EACH(vertex_a, shape_a.convex_vertices()) {
		dvec3 vertex_a_local = pos_b.to_local(pos_a.to_global(vertex_a));
		// TODO skip if current dist is smaller than distance between vertex_a and bounding sphere of B
		FOR_EACH(face_b, shape_b.faces())
            dist = std::min(dist, distance(face_b, vertex_a_local));
	}
	FOR_EACH(vertex_b, shape_b.convex_vertices()) {
		dvec3 vertex_b_local = pos_a.to_local(pos_b.to_global(vertex_b));
		// TODO skip if current dist is smaller than distance between vertex_b and bounding sphere of A
		FOR_EACH(face_a, shape_a.faces())
            dist = std::min(dist, distance(face_a, vertex_b_local));
	}

	// edge vs edge contacts
	FOR_EACH(edge_a, shape_a.convex_edges()) {
		segment3 edge_a_local = pos_b.to_local(pos_a.to_global(edge_a));
		// TODO skip if current dist is smaller than distance between edge_a and bounding sphere of B
		FOR_EACH(edge_b, shape_b.convex_edges())
            dist = std::min(dist, distance(edge_a_local, edge_b));
	}

    return dist;
}
#endif
