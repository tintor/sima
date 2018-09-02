#include "segment.h"
#include "aabb.h"
#include "exception.h"

// D - disjoint
// O - overlap
// V - vertex / vertex touch (could be colinear, but not overlapping)
// X - interior intersection
// A - T intersection: A or B is touching interior of PQ
// B - T intersection: P or Q is touching interior of AB

// TODO(marko) test that result will be stable if params are swapped (or if segments are reversed)
char relate(segment2 p, segment2 q) {
	int sp = Classify(p, q.a);
	int sq = Classify(p, q.b);
	if (sp == 0 && sq == 0) { // colinear
		if (Overlaps(aabb2(p), aabb2(q)))
			return 'O';
		if (Equals(p.a, q.a) || Equals(p.a, q.b) || Equals(p.b, q.a) || Equals(p.b, q.b))
			return 'V';
		return 'D';
	}

	if (Equals(p.a, q.a) || Equals(p.a, q.b) || Equals(p.b, q.a) || Equals(p.b, q.b))
		return 'V';

	int sab = Classify(q, p.a) * Classify(q, p.b);
	int spq = sp * sq;
	if (sab < 0 && spq < 0)
		return 'X';
	if (sab == 0 && spq < 0)
		return 'A';
	if (sab < 0 && spq == 0)
		return 'B';
	return 'D';
}

bool relate(segment2 p, ray2 q) {
	int sp = Classify(p, q.origin);
	int sq = Classify(p, q.origin + q.unit_dir);
	if (sp == 0 && sq == 0) { // colinear
		if (Overlaps(aabb2(p), aabb2(q.origin, q.infinity())))
			return 'O';
		if (Equals(p.a, q.origin) || Equals(p.b, q.origin))
			return 'V';
		return 'D';
	}

	if (Equals(p.a, q.origin) || Equals(p.b, q.origin))
		return 'V';

	int sab = Classify(q, p.a) * Classify(q, p.b);
	int spq = sp * sq;
	if (sab < 0 && spq < 0)
		return 'X';
	if (sab == 0 && spq < 0)
		return 'A';
	if (sab < 0 && spq == 0)
		return 'B';
	return 'D';
}

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

double distance(double4 a, double4 b) { return length(a - b); }

double distance(double4 a, segment3 b) { return distance(a, b.nearest(a)); }
double distance(segment3 a, double4 b) { return distance(b, a); }

double distance(segment3 a, segment3 b) {
    return sqrt(squared_distance(a, b));
}
