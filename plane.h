#pragma once
#include "format.h"

// TODO raw_plane: which stores non-unit normal

class plane {
public:
	static double4 init(double4 a, double4 b, double4 c) {
		double4 normal = normalize(compute_normal(a, b, c));
		return {normal.x, normal.y, normal.z, -dot(normal, a)};
	}

	plane() : n{NAN, NAN, NAN, NAN} { }
	plane(double4 normal, double d) : n{normal} { assert(is_unit(normal)); n.w = d; }
	plane(double4 a, double4 b, double4 c) : n{init(a, b, c)} { }
	plane(triangle3 v) : plane(v.a, v.b, v.c) { }

	plane(const plane& p) : n(p.n) { }
	plane& operator=(const plane& p) { n = p.n; return *this; }

	double4 normal() const { double4 e = n; e.w = 0; return e; }

	// Signed distance!
	double distance(double4 p) const {
		assert(p.w == 1.0);
		return dot(n, p);
	}

	// planes are 1um thick by default!
	int classify(double4 v, double eps = 0.5e-6) {
		double d = distance(v);
		if (d > eps) return 1;
		if (d < -eps) return -1;
		return 0;
	}

	// > 0, if P is on the positive side of plane ABC (right hand rule)
	static double sign(double4 a, double4 b, double4 c, double4 p) {
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

inline int Classify(plane p, double4 v) {
	return Classify(p.distance(v));
}
