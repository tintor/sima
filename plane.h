#pragma once
#include "format.h"

// TODO raw_plane: which stores non-unit normal

class plane {
public:
	static double4 init(double3 a, double3 b, double3 c) {
		double3 normal = normalize(compute_normal(a, b, c));
		return {normal.x, normal.y, normal.z, -dot(normal, a)};
	}

	plane() : n{NAN, NAN, NAN, NAN} { }
	plane(double3 normal, double d) : n{normal.x, normal.y, normal.z, d} { assert(is_unit(normal)); }
	plane(double3 a, double3 b, double3 c) : n{init(a, b, c)} { }
	plane(triangle3 v) : plane(v.a, v.b, v.c) { }

	plane(const plane& p) : n(p.n) { }
	plane& operator=(const plane& p) { n = p.n; return *this; }

	double3 normal() const { return n.xyz; }

	// Signed distance!
	double distance(double3 p) const {
		double4 p4 = d4(p);
		assert(p4.w == 1.0);
		return dot(n, p4);
	}

	// planes are 1um thick by default!
	int classify(double3 v, double eps = 0.5e-6) {
		double d = distance(v);
		if (d > eps) return 1;
		if (d < -eps) return -1;
		return 0;
	}

	// > 0, if P is on the positive side of plane ABC (right hand rule)
	static double sign(double3 a, double3 b, double3 c, double3 p) {
		return dot(cross(b - a, c - a), p - a);
	}

private:
	// n.xyz must be unit vector
	double4 n;
};

inline bool intersects(segment3 q, plane p, double eps) {
	double a = p.distance(q.a);
	double b = p.distance(q.b);
	return (a <= eps || b <= eps) && (a >= -eps || b >= -eps);
}

inline bool intersects(triangle3 q, plane p, double eps) {
	double a = p.distance(q.a);
	double b = p.distance(q.b);
	double c = p.distance(q.c);
	return (a <= eps || b <= eps || c <= eps) && (a >= -eps || b >= -eps || c >= -eps);
}
