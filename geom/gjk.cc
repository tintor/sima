#include <core/exception.h>
#include <geom/gjk.h>
#include <geom/segment.h>

constexpr double2 ZERO2 = {0, 0};
constexpr double4 ZERO3 = {0, 0, 0, 0};

// are P and Q of different sides of line AB?
inline bool on_different_sides(double2 a, double2 b, double2 p, double2 q) {
	double sp = signed_double_area(a, b, p);
	double sq = signed_double_area(a, b, q);
	return sp * sq < 0;
}

bool gjk_simplex(static_vector<double2, 3>& simplex, double2& dir) {
	if (simplex.size() == 1) {
		double2 a = simplex[0];
		dir = -a;
		return false;
	}
	if (simplex.size() == 2) {
		double2 b = simplex[0], a = simplex[1];
		double2 ab = b - a;
		if (dot(ab, -a) > 0) {
			double2 n = {ab.y, -ab.x};
			if (dot(n, -a) > 0)
				dir = n;
			else
				dir = -n;
		} else {
			dir = -a;
		}
		return false;
	}
	if (simplex.size() == 3) {
		double2 c = simplex[0], b = simplex[1], a = simplex[2];
		if (on_different_sides(a, b, c, ZERO2)) {
			simplex = {b, a};
			return gjk_simplex(simplex, dir);
		}
		if (on_different_sides(a, c, b, ZERO2)) {
			simplex = {c, a};
			return gjk_simplex(simplex, dir);
		}
		return true;
	}
	THROW(runtime_error);
}

// are P and Q on different sides of plane ABC?
inline bool on_different_sides(double4 a, double4 b, double4 c, double4 p, double4 q) {
	double4 n = cross(b - a, c - a);
	return dot(n, p - a) * dot(n, q - a) < 0;
}

bool gjk_simplex(static_vector<double4, 4>& simplex, double4& dir) {
	if (simplex.size() == 1) {
		double4 a = simplex[0];
		dir = -a;
		return false;
	}
	if (simplex.size() == 2) {
		double4 b = simplex[0], a = simplex[1];
		double4 ab = b - a;
		dir = (dot(ab, -a) > 0) ? cross(cross(ab, -a), ab) : -a;
		return false;
	}
	if (simplex.size() == 3) {
		double4 c = simplex[0], b = simplex[1], a = simplex[2];
		double4 up = a + cross(c - a, b - a);
		if (on_different_sides(a, b, up, c, ZERO3)) {
			simplex = {b, a};
			return gjk_simplex(simplex, dir);
		}
		if (on_different_sides(a, c, up, b, ZERO3)) {
			simplex = {c, a};
			return gjk_simplex(simplex, dir);
		}
		dir = (dot(up, -a) > 0) ? up : -up;
		return false;
	}
	if (simplex.size() == 4) {
		double4 d = simplex[0], c = simplex[1], b = simplex[2], a = simplex[3];
		if (on_different_sides(a, b, c, d, ZERO3)) {
			simplex = {c, b, a};
			return gjk_simplex(simplex, dir);
		}
		if (on_different_sides(a, b, d, c, ZERO3)) {
			simplex = {d, b, a};
			return gjk_simplex(simplex, dir);
		}
		if (on_different_sides(a, c, d, b, ZERO3)) {
			simplex = {d, c, a};
			return gjk_simplex(simplex, dir);
		}
		return true;
	}
	THROW(runtime_error);
}
