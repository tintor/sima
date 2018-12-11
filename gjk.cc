#include "gjk.h"
#include "exception.h"

struct Transform3 {
};

class Shape3 {
	enum class Kind { Sphere, ConeSphere, Capsule, Polyhedron, Box, RoundedBox } kind;
	double r1, r2;
	double4 size;
	vector<double4> vertices; // TODO reduce space overlap vector<> with

	double4 support(double4 dir);
};

/*
double ß = 0;
double Θ = 2;
double Ω = 3;
double Ψ = 1;
double Δ = 0;
*/

double4 Shape3::support(double4 dir) {

	switch (kind) {
	case Kind::Sphere:
		return dir * (r1 / length(dir));

	case Kind::ConeSphere:
		if (dir.x > 0)
			return double4{size.x, 0, 0, 0} + dir * (r1 / length(dir));
		return double4{-size.x, 0, 0, 0};

	case Kind::Capsule:
		if (dir.x > 0)
			return double4{size.x, 0, 0, 0} + dir * (r1 / length(dir));
		return double4{-size.x, 0, 0, 0} + dir * (r2 / length(dir));

	case Kind::Polyhedron: {
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

// TODO generalize support function to make it work for
// mix of point, box, cylinder-capsule, cone-capsule, sphere, convex polyhedron, rounded convex polyhedron
// capsule support is easy: convex hull of two spheres!
// TODO support any number of spheres with radius per sphere!
// TODO make it work with untransformed shapes!
// TODO support of box can be computed faster than for 8 vertex convex polyhedron

static double4 SphereSupport(double4 center, double radius, double4 dir) {
	return center + dir * (radius / length(dir));
}

template<typename V>
static V Support(span<const V> vertices, V dir) {
	double md = -std::numeric_limits<double>::infinity();
	V mv;
	for (V v : vertices) {
		double d = dot(v, dir);
		if (d > md) {
			md = d;
			mv = v;
		}
	}
	return mv;
}

// ====
//  2D
// ====

static double2 NearestSimplex2(double2 b, double2 a) {
	double2 ab = b - a;
	if (dot(ab, -a) > 0)
		return cross(cross(ab, -a), ab);
	return -a;
}

// returns true if contains origin
static bool NearestSimplex3(double2 c, double2 b, double2 a) {
	THROW(not_implemented);
}

// TODO if no axis then difference of centers is better guess!
// TODO it is possible to terminate early for Intersects
// TODO return final axis
bool AreConvexHullsIntersecting(span<const double2> p, span<const double2> q, double2* axis) {
	double2 initial_axis = axis ? *axis : (q[0] - p[0]);
	double2 A = Support(p, initial_axis) - Support(q, -initial_axis);
	double2 D = -A;
	array<double2, 3> simplex = {A};

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
	return NearestSimplex3(simplex[0], simplex[1], simplex[2]);
}

double DistanceBetweenConvexHulls(span<const double2> p, span<const double2> q, double2* axis) {
	THROW(not_implemented);
}

double SignedDistanceBetweenConvexHulls(span<const double2> p, span<const double2> q, double2* axis) {
	THROW(not_implemented);
}

optional<segment2> MinSegmentBetweenConvexHulls(span<const double2> p, span<const double2> q, double2* axis) {
	THROW(not_implemented);
}

// ====
//  3D
// ====

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
bool AreConvexHullsIntersecting(span<const double4> p, span<const double4> q, double4* axis) {
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
}

double DistanceBetweenConvexHulls(span<const double4> p, span<const double4> q, double4* axis) {
	THROW(not_implemented);
}

double SignedDistanceBetweenConvexHulls(span<const double4> p, span<const double4> q, double4* axis) {
	THROW(not_implemented);
}

optional<segment3> MinSegmentBetweenConvexHulls(span<const double4> p, span<const double4> q, double4* axis) {
	THROW(not_implemented);
}
