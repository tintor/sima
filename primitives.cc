#include "primitives.h"

#ifdef xxx
// infinite line vs. point
real line_point_squared_distance(dvec3 a, dvec3 b, dvec3 p) {
	dvec3 ba = b - a, pa = p - a;
    real t = dot(pa, ba) / dot(ba, ba);
    return squared(pa - ba * t);
}

// line segment
std::pair<segment3, segment3::NearestCase> segment3::Nearest(segment3 p, segment3 q) {
	dvec3 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	real aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr real inf = std::numeric_limits<real>::max();
	constexpr real tiny = 1e-8;

	// ray/ray
	real d = aa * bb - ab * ab;
	real s = ab * bc - bb * ac;
	real t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0))
		return std::pair<segment3, NearestCase>(segment3(p.a + A * (s / d), q.a + B * (t / d)), NearestCase::RayRay);

	// ray/endpoint
	real s0 = (aa >= tiny) ? -ac / aa : -1;
	real s1 = (aa >= tiny) ? (ab - ac) / aa : -1;
	real t0 = (bb >= tiny) ? bc / bb : -1;
	real t1 = (bb >= tiny) ? (ab + bc) / bb : -1;

	real d1 = (0 <= s0 && s0 <= 1) ? squared(C + A*s0) : inf;
	real d2 = (0 <= s1 && s1 <= 1) ? squared(C - B + A*s1) : inf;
	real d3 = (0 <= t0 && t0 <= 1) ? squared(B*t0 - C) : inf;
	real d4 = (0 <= t1 && t1 <= 1) ? squared(B*t1 - C - A) : inf;

	// endpoint/endpoint
	real d5 = squared(C);
	real d6 = squared(C + A);
	real d7 = squared(C - B);
	real d8 = squared(C + A - B);

	real dm = std::min(min(d1, d2, d3, d4), min(d5, d6, d7, d8));

	if (d1 == dm)
		return std::pair(segment3(p.a + A*s0, q.a), NearestCase::RayPoint);
	if (d2 == dm)
		return std::pair(segment3(p.a + A*s1, q.b), NearestCase::RayPoint);
	if (d3 == dm)
		return std::pair(segment3(p.a, q.a + B*t0), NearestCase::RayPoint);
	if (d4 == dm)
		return std::pair(segment3(p.b, q.a + B*t1), NearestCase::RayPoint);
	if (d5 == dm)
		return std::pair(segment3(p.a, q.a), NearestCase::PointPoint);
	if (d6 == dm)
		return std::pair(segment3(p.b, q.a), NearestCase::PointPoint);
	if (d7 == dm)
		return std::pair(segment3(p.a, q.b), NearestCase::PointPoint);
	if (d8 == dm)
		return std::pair(segment3(p.b, q.b), NearestCase::PointPoint);

	throw new std::exception();
}

real segment3::squared_distance(segment3 p, segment3 q) {
	dvec3 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	real aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr real inf = std::numeric_limits<real>::max();
	constexpr real tiny = 1e-8;

	// ray/ray
	real d = aa * bb - ab * ab;
	real s = ab * bc - bb * ac;
	real t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0))
		return squared(C + A * (s / d) - B * (t / d));

	// ray/endpoint
	real s0 = (aa >= tiny) ? -ac / aa : -1;
	real s1 = (aa >= tiny) ? (ab - ac) / aa : -1;
	real t0 = (bb >= tiny) ? bc / bb : -1;
	real t1 = (bb >= tiny) ? (ab + bc) / bb : -1;

	real d1 = (0 <= s0 && s0 <= 1) ? squared(C + A*s0) : inf;
	real d2 = (0 <= s1 && s1 <= 1) ? squared(C - B + A*s1) : inf;
	real d3 = (0 <= t0 && t0 <= 1) ? squared(B*t0 - C) : inf;
	real d4 = (0 <= t1 && t1 <= 1) ? squared(B*t1 - C - A) : inf;

	// endpoint/endpoint
	real d5 = squared(C);
	real d6 = squared(C + A);
	real d7 = squared(C - B);
	real d8 = squared(C + A - B);

	return std::min(min(d1, d2, d3, d4), min(d5, d6, d7, d8));
}

real distance(const dvec3& a, const dvec3& b) { return l2Norm(a - b); }

real distance(const dvec3& a, const segment3& b) { return distance(a, b.nearest(a)); }
real distance(const segment3& a, const dvec3& b) { return distance(b, a); }

real distance(const segment3& a, const segment3& b) {
    return sqrt(segment3::squared_distance(a, b));
}

constexpr bool inside_triangle_prism(const dvec3& p, const triangle3& m, const dvec3& normal) {
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + normal, p) > 0)
			return false;
	return true;
}

real distance(const dvec3& p, const triangle3& m) {
    dvec3 n = Normal(m);
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + n, p) > 0)
			return distance(p, segment3(*a, *b));
	return abs(dot(n, p - m.a));
}

real distance(const dvec3& v, const triangle3& m, const plane& p) {
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + p.normal, v) > 0)
			return distance(v, segment3(*a, *b));
	return abs(p.distance(v));
}

real distance(const triangle3& m, const dvec3& v) { return distance(v, m); }

// from RealTimeCollisionDetection book
bool intersects(const line3& e, const triangle3& m) {
	auto d = e.b - e.a;
	auto pa = m.a - e.a;
	auto pb = m.b - e.a;
	auto pc = m.c - e.a;

	auto n = cross(d, pc);
	return dot(pb, n) >= 0 && dot(pa, n) /*intentional*/<= 0 && dot(cross(d, pb), pa) >= 0;
}

// from RealTimeCollisionDetection book
bool intersects2(const line3& e, const triangle3& m) {
	auto d = e.b - e.a;
	auto n = cross(d, e.b);
	auto s = dot(n, m.c - m.b);
	auto t = dot(n, m.a - m.c);
	// TODO cross products can be precomputed and triangle replaced with plucker!
	return dot(d, cross(m.c, m.b)) + s >= 0 && dot(d, cross(m.a, m.c)) + t >= 0 && dot(d, cross(m.b, m.a)) - s - t >= 0;
}

// from RealTimeCollisionDetection book
bool intersection(const line3& e, const triangle3& m, /*out*/vec3& result) {
	auto d = e.b - e.a;
	auto pa = m.a - e.a;
	auto pb = m.b - e.a;
	auto pc = m.c - e.a;

	auto n = cross(d, pc);
	auto u = dot(pb, n);
	if (u < 0)
		return false;
	auto v = -dot(pa, n);
	if (v < 0)
	 	return false;
	auto w = dot(cross(d, pb), pa);
	if (w < 0)
		return false;

	auto denom = u + v + w;
	if (denom < 1e-8)
		throw new std::runtime_error("planar case");
	result = (m.a * u + m.b * v + m.c * w) / denom;
	return true;
}

// Note: Ignores coplanar case!
bool intersects_in_point(const segment3& e, const triangle3& m) {
	auto d = e.b - e.a;
	auto pa = m.a - e.a;
	auto pb = m.b - e.a;
	auto pc = m.c - e.a;

	auto n = cross(d, pc);
	auto u = dot(pb, n);
	if (u < 0)
		return false;
	auto v = -dot(pa, n);
	if (v < 0)
	 	return false;
	auto w = dot(cross(d, pb), pa);
	if (w < 0)
		return false;

	// Note: rejecting parallel case
	auto denom = u + v + w;
	if (denom < 1e-8)
		return false;

	// G is intersection of line with triangle
	auto g = (m.a * u + m.b * v + m.c * w) / denom;

	// Check if G is on the segment
	auto t = dot(g - e.a, d);
	return 0 <= t && t <= dot(d, d);
}

// Note: Ignores coplanar case!
bool intersects_in_point(const ray3& e, const triangle3& m) {
	auto pa = m.a - e.origin;
	auto pb = m.b - e.origin;
	auto pc = m.c - e.origin;

	auto n = cross(e.dir, pc);
	auto u = dot(pb, n);
	if (u < 0)
		return false;
	auto v = -dot(pa, n);
	if (v < 0)
	 	return false;
	auto w = dot(cross(e.dir, pb), pa);
	if (w < 0)
		return false;

	// Note: rejecting parallel case
	auto denom = u + v + w;
	if (denom < 1e-8)
		return false;

	// G is intersection of line with triangle
	auto g = (m.a * u + m.b * v + m.c * w) / denom;

	// Check if G is on the ray
	return 0 <= dot(g - e.origin, e.dir);
}

real disjoint_distance(const segment3& e, const triangle3& m) {
	// Needed to handle a single point intersection case.
	// Not needed for case when segment overlaps with triangle (coplanar case).
	if (intersects_in_point(e, m))
		return 0.0;

    auto d1 = segment3::squared_distance(e, segment3(m.a, m.b));
    auto d2 = segment3::squared_distance(e, segment3(m.b, m.c));
    auto d3 = segment3::squared_distance(e, segment3(m.c, m.a));
    auto d = sqrt(min(d1, d2, d3));

    auto n = Normal(m);
    if (inside_triangle_prism(e.a, m, n))
        d = std::min(d, abs(dot(n, e.a - m.a)));
    if (inside_triangle_prism(e.b, m, n))
        d = std::min(d, abs(dot(n, e.b - m.a)));
    return d;
}

real distance(const segment3& e, const triangle3& m) {
	// Needed to handle a single point intersection case.
	// Not needed for case when segment overlaps with triangle (coplanar case).
	if (intersects_in_point(e, m))
		return 0.0;

    return disjoint_distance(e, m);
}

real distance(const triangle3 m, const segment3& e) { return distance(e, m); }

// Assuming these two objects do not intersect!
real disjoint_distance(const triangle3& p, const triangle3& q) {
    real d = std::numeric_limits<real>::max();
    // TODO only compute these 6 distances if inside triangle prisms
    FOR(i, 3)
        d = min(d, distance(p, q[i]), distance(p[i], q));

	real ds = 1e100;
	FOR_EACH_EDGE(pa, pb, p)
		FOR_EACH_EDGE(qa, qb, q)
			ds = std::min(ds, segment3::squared_distance(segment3(*pa, *pb), segment3(*qa, *qb)));
    return std::min(d, sqrt(ds));
}

real distance(const triangle3& p, const triangle3& q) {
	// Needed to handle a single point intersection case.
	// Not needed for case when one triangle overlaps with another (coplanar case).
	// TODO here it is worthwhile to precompute cross products of Q
	FOR_EACH_EDGE(pa, pb, p)
		if (intersects_in_point(segment3(*pa, *pb), q))
			return 0.0;

	return disjoint_distance(p, q);
}

// Angle between oriented triangles ABC and BAD.
// If edge is convex then angle will be <PI
// If edge is planar then angle will be =PI
// If edge is concave than angle will be >PI
real edge_angle(const dvec3& a, const dvec3& b, const dvec3& c, const dvec3& d) {
	throw new std::runtime_error("edge_angle() unimplemented");
}

sphere merge_spheres(const sphere& a, const sphere& b) {
	auto d = b.center - a.center;
	auto d2 = squared(d);

	if (squared(a.radius - b.radius) >= d2)
		return (a.radius > b.radius) ? a : b;

	auto dist = sqrt(d2);
	sphere c;
	c.radius = (a.radius + b.radius + dist) / 2;
	// division is safe here for small d, as center will converge to a.center
	c.center = a.center + ((c.radius - a.radius) / dist) * d;
	return c;
}
#endif
