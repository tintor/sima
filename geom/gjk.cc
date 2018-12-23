#include <geom/gjk.h>
#include <geom/buffer.h>
#include <geom/pose.h>
#include <core/range.h>

// double ß = 0;
// double Θ = 2;
// double Ω = 3;
// double Ψ = 1;
// double Δ = 0;

// ====
//  2D
// ====

/*class Shape2 {
	enum class Kind { Circle, Capsule, Polygon, Box, RoundedBox } kind;
	double r1, r2;
	double2 size;
	vector<double2> vertices; // TODO reduce space overlap vector<> with

	double2 support(double2 dir);
};

double2 Shape2::support(double2 dir) {
	switch (kind) {
	case Kind::Circle:
		return dir * (r1 / length(dir));

	case Kind::Capsule:
		if (dir.x > 0)
			return double2{size.x, 0} + dir * (r1 / length(dir));
		return double2{-size.x, 0} + dir * (r2 / length(dir));

	case Kind::Polygon: {
		// TODO could also convert dir to angle and do a table lookup
		double md = -std::numeric_limits<double>::infinity();
		double2 mv;
		for (double2 v : vertices) {
			double d = dot(v, dir);
			if (d > md) {
				md = d;
				mv = v;
			}
		}
		if (r1 != 0)
			mv += dir * (r1 / length(dir));
		return mv;
	}
	case Kind::Box:
		return sign_no_zero(dir) * size;

	case Kind::RoundedBox:
		return sign_no_zero(dir) * size + dir * (r1 / length(dir));
	}
	THROW(runtime_error);
}*/

static double2 NearestSimplex2(double2 b, double2 a) {
	double2 ab = b - a;
	if (dot(ab, -a) > 0)
		return cross(cross(ab, -a), ab);
	return -a;
}

// returns true if contains origin
static bool NearestSimplex3(double2 c, double2 b, double2 a) {
	double sz = signed_double_area(a, b, 0);
	double sc = signed_double_area(a, b, c);
	if (sz * sc < 0)
		return false;

	sz = signed_double_area(a, c, 0);
	double sb = signed_double_area(a, c, b);
	return sz * sb >= 0;
}

// center of mass is assumed to be at zero
struct Shape2 {
	double radius;
	vector<double2> vertices;
};

pair<double, double> MinMaxDotProduct(cspan<double2> p, double2 d) {
	double v_max = -std::numeric_limits<double>::infinity();
	double v_min = std::numeric_limits<double>::infinity();
	for (double2 v : p) {
		double vd = dot(v, d);
		maximize(v_max, vd);
		minimize(v_min, vd);
	}
	return {v_min, v_max};
}

// signed distance between two intervals
double Distance(pair<double, double> a, pair<double, double> b) {
	return max(b.first - a.second, a.first - b.second);
}

// mid point between two intervals
double Midpoint(pair<double, double> a, pair<double, double> b) {
	if (b.first - a.second > a.first - b.second)
		return (b.first + a.second) / 2;
	return (a.second + b.first) / 2;
}

double2 CentroidOfConvexIntersection(const Shape2& sa, Pose2 pa, const Shape2& sb, Pose2 pb) {
	THROW(not_implemented);
}

pair<double, double> IntersectionInterval(const Shape2& s, double2 dir, double2 p) {
	THROW(not_implemented);
}

struct Manifold2 {
	double depth;
	double2 normal; // in world coordinates pointing towards A
	double2 centroid; // only if intersecting
};

// robot arm demo:
// min: box vs box
// better: convex_poly vs convex_poly

/*optional<Manifold2> AreConvexShapesIntersecting(Pose2 pa, const Shape2& sa, Pose2 pb, const Shape2& sb) {
	double distance = -INF;
	double2 dir2;
	constexpr int N = 120;
	for (int i : range(N)) {
		double e = (PI / N) * i;
		double2 dir{cos(e), sin(e)}; // TODO precompute
		auto ia = MinMaxDotProduct(sa.vertices, dir);
		auto ib = MinMaxDotProduct(sb.vertices, dir);
		double d = Distance(ia, ib);
		if (d > 0)
			return nullopt;
		if (d > distance) {
			distance = d;
			dir2 = dir;
		}
	}

	// TODO intersect boundaries
	// if N == 2 return normal based on two points
	// if N == 4 return normal based on dir
		//
	return comp;
}*/

// TODO if no axis then difference of centers is better guess!
// TODO it is possible to terminate early for Intersects
// TODO return final axis
/*bool AreConvexHullsIntersecting(cspan<double2> p, cspan<double2> q, double2* axis) {
	double2 D = axis ? *axis : (q[0] - p[0]);
	double2 A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;
		return false;
	}

	double2 S1 = A;
	D = -A;
	A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;
		return false;
	}

	double2 S2 = A;
	D = NearestSimplex2(S1, S2);
	A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;
		return false;
	}

	return NearestSimplex3(S1, S2, A);
}*/

/*double DistanceBetweenConvexHulls(cspan<double2> p, cspan<double2> q, double2* axis) {
	double2 initial_axis = axis ? *axis : (q[0] - p[0]);
	double2 A = Support(p, initial_axis) - Support(q, -initial_axis);
	double2 D = -A;

	double2 S1 = A;
	A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;

		// definitely not intersecting
		D = -A;
		while (true) {
			A = Support(p, D) - Support(q, -D);
			if (A == -D)
				return length(A);
			D = -A;
		}
	}

	double2 S2 = A;
	D = NearestSimplex2(S1, S2);

	A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;

		// definitely not intersecting
		D = -A;
		while (true) {
			A = Support(p, D) - Support(q, -D);
			if (A == -D)
				return length(A);
			D = -A;
		}
	}

	if (NearestSimplex3(S1, S2, A))
		return 0;

	// definitely not intersecting
	D = -A; // TODO set D from Simplex3
	while (true) {
		A = Support(p, D) - Support(q, -D);
		if (A == -D)
			return length(A);
		D = -A;
	}
}*/

double SignedDistanceBetweenConvexHulls(cspan<double2> p, cspan<double2> q, double2* axis) {
	THROW(not_implemented);
}

optional<segment2> MinSegmentBetweenConvexHulls(cspan<double2> p, cspan<double2> q, double2* axis) {
	THROW(not_implemented);
}

/*double DistanceBetweenConvexHulls2(cspan<double2> p, cspan<double2> q, int pi, int qi) {
	// start with two candidate points
	// try to jiggle (move each point prev/next as long as it is reducing distance)
	// when we have two closest points, find which of the 4 vertex/edge is closest

	// TODO how to detect intersection?
}*/

// ====
//  3D
// ====

class Shape3 {
	enum class Kind { Sphere, Capsule, Polyhedron, Box, RoundedBox } kind;
	double r1, r2;
	double4 size;
	vector<double4> vertices; // TODO reduce space overlap vector<> with

	double4 support(double4 dir);
};

double4 Shape3::support(double4 dir) {
	switch (kind) {
	case Kind::Sphere:
		return dir * (r1 / length(dir));

	case Kind::Capsule:
		if (dir.x > 0)
			return double4{size.x, 0, 0, 1} + dir * (r1 / length(dir));
		return double4{-size.x, 0, 0, 1} + dir * (r2 / length(dir));

	case Kind::Polyhedron: {
		// TODO could also construct a decision tree for large number of vertices!
		double md = -std::numeric_limits<double>::infinity();
		double4 mv;
		for (double4 v : vertices) {
			double d = dot(v, dir);
			if (d > md) {
				md = d;
				mv = v;
			}
		}
		if (r1 != 0)
			mv += dir * (r1 / length(dir));
		return mv;
	}
	case Kind::Box:
		return sign_no_zero(dir) * size;

	case Kind::RoundedBox:
		return sign_no_zero(dir) * size + dir * (r1 / length(dir));
	}
	THROW(runtime_error);
}

int SignedDistance(const Shape3& p, const Shape3& q, double4& axis) {
	THROW(not_implemented);
}

static double4 NearestSimplex2(double4 b, double4 a) {
	double4 ab = b - a;
	double4 ao = -a;

	if (dot(ab, ao) > 0)
		return cross(cross(ab, ao), ab);
	return ao;
}

static double4 NearestSimplex3(double4 c, double4 b, double4 a) {
	THROW(not_implemented);
}

// returns true if contains origin
static bool NearestSimplex4(double4 d, double4 c, double4 b, double4 a) {
	THROW(not_implemented);
}

// TODO if no axis then difference of centers is better guess!
// TODO it is possible to terminate early for Intersects
// TODO return final axis
/*bool AreConvexHullsIntersecting(cspan<double4> p, cspan<double4> q, double4* axis) {
	double4 initial_axis = axis ? *axis : (q[0] - p[0]);
	double4 A = Support(p, initial_axis) - Support(q, -initial_axis);
	double4 D = -A;
	array<double4, 4> simplex = {A};

	A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;
		return false;
	}

	simplex[1] = A;
	D = NearestSimplex2(simplex[0], simplex[1]);

	A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;
		return false;
	}

	simplex[2] = A;
	D = NearestSimplex3(simplex[0], simplex[1], simplex[2]);

	A = Support(p, D) - Support(q, -D);
	if (dot(A, D) < 0) {
		if (axis)
			*axis = D;
		return false;
	}

	simplex[3] = A;
	return NearestSimplex4(simplex[0], simplex[1], simplex[2], simplex[3]);
}*/

double DistanceBetweenConvexHulls(cspan<double4> p, cspan<double4> q, double4* axis) {
	THROW(not_implemented);
}

double SignedDistanceBetweenConvexHulls(cspan<double4> p, cspan<double4> q, double4* axis) {
	THROW(not_implemented);
}

optional<segment3> MinSegmentBetweenConvexHulls(cspan<double4> p, cspan<double4> q, double4* axis) {
	THROW(not_implemented);
}
