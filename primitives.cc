#include "primitives.h"
#include "exception.h"

// line segment
pair<segment3, NearestCase> nearest(segment3 p, segment3 q) {
	double4 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	double aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr double inf = std::numeric_limits<double>::max();
	constexpr double tiny = 1e-8;

	// ray/ray
	double d = aa * bb - ab * ab;
	double s = ab * bc - bb * ac;
	double t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0))
		return pair<segment3, NearestCase>(segment3(p.a + A * (s / d), q.a + B * (t / d)), NearestCase::RayRay);

	// ray/endpoint
	double s0 = (aa >= tiny) ? -ac / aa : -1;
	double s1 = (aa >= tiny) ? (ab - ac) / aa : -1;
	double t0 = (bb >= tiny) ? bc / bb : -1;
	double t1 = (bb >= tiny) ? (ab + bc) / bb : -1;

	double d1 = (0 <= s0 && s0 <= 1) ? squared(C + A*s0) : inf;
	double d2 = (0 <= s1 && s1 <= 1) ? squared(C - B + A*s1) : inf;
	double d3 = (0 <= t0 && t0 <= 1) ? squared(B*t0 - C) : inf;
	double d4 = (0 <= t1 && t1 <= 1) ? squared(B*t1 - C - A) : inf;

	// endpoint/endpoint
	double d5 = squared(C);
	double d6 = squared(C + A);
	double d7 = squared(C - B);
	double d8 = squared(C + A - B);

	double dm = std::min(min(d1, d2, d3, d4), min(d5, d6, d7, d8));

	if (d1 == dm)
		return pair(segment3(p.a + A*s0, q.a), NearestCase::RayPoint);
	if (d2 == dm)
		return pair(segment3(p.a + A*s1, q.b), NearestCase::RayPoint);
	if (d3 == dm)
		return pair(segment3(p.a, q.a + B*t0), NearestCase::RayPoint);
	if (d4 == dm)
		return pair(segment3(p.b, q.a + B*t1), NearestCase::RayPoint);
	if (d5 == dm)
		return pair(segment3(p.a, q.a), NearestCase::PointPoint);
	if (d6 == dm)
		return pair(segment3(p.b, q.a), NearestCase::PointPoint);
	if (d7 == dm)
		return pair(segment3(p.a, q.b), NearestCase::PointPoint);
	if (d8 == dm)
		return pair(segment3(p.b, q.b), NearestCase::PointPoint);

	THROW(runtime_error, "unreachable");
}

double squared_distance(segment3 p, segment3 q) {
	double4 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	double aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr double inf = std::numeric_limits<double>::max();
	constexpr double tiny = 1e-8;

	// ray/ray
	double d = aa * bb - ab * ab;
	double s = ab * bc - bb * ac;
	double t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0))
		return squared(C + A * (s / d) - B * (t / d));

	// ray/endpoint
	double s0 = (aa >= tiny) ? -ac / aa : -1;
	double s1 = (aa >= tiny) ? (ab - ac) / aa : -1;
	double t0 = (bb >= tiny) ? bc / bb : -1;
	double t1 = (bb >= tiny) ? (ab + bc) / bb : -1;

	double d1 = (0 <= s0 && s0 <= 1) ? squared(C + A*s0) : inf;
	double d2 = (0 <= s1 && s1 <= 1) ? squared(C - B + A*s1) : inf;
	double d3 = (0 <= t0 && t0 <= 1) ? squared(B*t0 - C) : inf;
	double d4 = (0 <= t1 && t1 <= 1) ? squared(B*t1 - C - A) : inf;

	// endpoint/endpoint
	double d5 = squared(C);
	double d6 = squared(C + A);
	double d7 = squared(C - B);
	double d8 = squared(C + A - B);

	return min(d1, d2, d3, d4, d5, d6, d7, d8);
}

double squared_distance(line3 e, double4 p) {
	double4 pa = p - e.origin;
    return squared(pa - e.unit_dir * dot(pa, e.unit_dir));
}

// TODO
double squared_distance(segment3 p, segment3 q);
double squared_distance(line3 m, double4 p);

double distance(double4 a, double4 b) { return length(a - b); }

double distance(double4 a, segment3 b) { return distance(a, b.nearest(a)); }
double distance(segment3 a, double4 b) { return distance(b, a); }

double distance(segment3 a, segment3 b) {
    return sqrt(squared_distance(a, b));
}

constexpr bool inside_triangle_prism(double4 p, triangle3 m, double4 normal) {
	for (auto [a, b] : edgesOf(m))
		if (plane::sign(a, b, a + normal, p) > 0)
			return false;
	return true;
}

double distance(double4 p, triangle3 m) {
    double4 normal = compute_normal(m);
	for (auto [a, b] : edgesOf(m))
		if (plane::sign(a, b, a + normal, p) > 0)
			return distance(p, segment3(a, b));
	return abs(dot(normal, p - m.a));
}

double distance(double4 v, triangle3 m, const plane& p) {
	for (auto [a, b] : edgesOf(m))
		if (plane::sign(a, b, a + p.normal(), v) > 0)
			return distance(v, segment3(a, b));
	return abs(p.distance(v));
}

double distance(triangle3 m, double4 v) { return distance(v, m); }

// from RealTimeCollisionDetection book
bool intersects(line3 e, triangle3 m) {
	auto pa = m.a - e.origin;
	auto pb = m.b - e.origin;
	auto pc = m.c - e.origin;

	auto n = cross(e.unit_dir, pc);
	return dot(pb, n) >= 0 && dot(pa, n) /*intentional*/<= 0 && dot(cross(e.unit_dir, pb), pa) >= 0;
}

// from RealTimeCollisionDetection book
bool intersects2(line3 e, triangle3 m) {
	auto d = e.unit_dir;
	auto n = cross(d, e.origin + d);
	auto s = dot(n, m.c - m.b);
	auto t = dot(n, m.a - m.c);
	// TODO cross products can be precomputed and triangle replaced with plucker!
	return dot(d, cross(m.c, m.b)) + s >= 0 && dot(d, cross(m.a, m.c)) + t >= 0 && dot(d, cross(m.b, m.a)) - s - t >= 0;
}

// from RealTimeCollisionDetection book
bool intersection(line3 e, triangle3 m, /*out*/double4& result) {
	auto pa = m.a - e.origin;
	auto pb = m.b - e.origin;
	auto pc = m.c - e.origin;

	auto n = cross(e.unit_dir, pc);
	auto u = dot(pb, n);
	if (u < 0)
		return false;
	auto v = -dot(pa, n);
	if (v < 0)
	 	return false;
	auto w = dot(cross(e.unit_dir, pb), pa);
	if (w < 0)
		return false;

	auto denom = u + v + w;
	if (denom < 1e-8)
		THROW(runtime_error, "planar case");
	result = (m.a * u + m.b * v + m.c * w) / denom;
	return true;
}

// Note: Ignores coplanar case!
bool intersects_in_point(segment3 e, triangle3 m) {
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
bool intersects_in_point(ray3 e, triangle3 m) {
	auto pa = m.a - e.origin;
	auto pb = m.b - e.origin;
	auto pc = m.c - e.origin;

	auto n = cross(e.unit_dir, pc);
	auto u = dot(pb, n);
	if (u < 0)
		return false;
	auto v = -dot(pa, n);
	if (v < 0)
	 	return false;
	auto w = dot(cross(e.unit_dir, pb), pa);
	if (w < 0)
		return false;

	// Note: rejecting parallel case
	auto denom = u + v + w;
	if (denom < 1e-8)
		return false;

	// G is intersection of line with triangle
	auto g = (m.a * u + m.b * v + m.c * w) / denom;

	// Check if G is on the ray
	return 0 <= dot(g - e.origin, e.unit_dir);
}

double disjoint_distance(segment3 e, triangle3 m) {
	// Needed to handle a single point intersection case.
	// Not needed for case when segment overlaps with triangle (coplanar case).
	if (intersects_in_point(e, m))
		return 0.0;

    auto d1 = squared_distance(e, segment3(m.a, m.b));
    auto d2 = squared_distance(e, segment3(m.b, m.c));
    auto d3 = squared_distance(e, segment3(m.c, m.a));
    auto d = sqrt(min(d1, d2, d3));

    auto n = compute_normal(m);
    if (inside_triangle_prism(e.a, m, n))
        d = min(d, abs(dot(n, e.a - m.a)));
    if (inside_triangle_prism(e.b, m, n))
        d = min(d, abs(dot(n, e.b - m.a)));
    return d;
}

double distance(segment3 e, triangle3 m) {
	// Needed to handle a single point intersection case.
	// Not needed for case when segment overlaps with triangle (coplanar case).
	if (intersects_in_point(e, m))
		return 0.0;

    return disjoint_distance(e, m);
}

double distance(const triangle3 m, segment3 e) { return distance(e, m); }

// Assuming these two objects do not intersect!
double disjoint_distance(triangle3 p, triangle3 q) {
    double d = std::numeric_limits<double>::max();
    // TODO only compute these 6 distances if inside triangle prisms
    for (auto i : range(3))
        d = min(d, distance(p, q[i]), distance(p[i], q));

	double ds = 1e100;
	for (auto ep : edgesOf(p))
		for (auto eq : edgesOf(q))
			ds = min(ds, squared_distance(ep, eq));
    return min(d, sqrt(ds));
}

double distance(triangle3 p, triangle3 q) {
	// Needed to handle a single point intersection case.
	// Not needed for case when one triangle overlaps with another (coplanar case).
	// TODO here it is worthwhile to precompute cross products of Q
	for (auto e : edgesOf(p))
		if (intersects_in_point(e, q))
			return 0.0;

	return disjoint_distance(p, q);
}

// Angle between oriented triangles ABC and BAD.
// If edge is convex then angle will be <PI
// If edge is planar then angle will be =PI
// If edge is concave than angle will be >PI
double edge_angle(double4 a, double4 b, double4 c, double4 d) {
	THROW(not_implemented);
}
