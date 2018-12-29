#include <core/std.h>
#include <core/range.h>
#include <core/format.h>
#include <core/dynamic_array.h>
#include <geom/vector.h>

bool gjk_simplex(static_vector<double2, 3>& simplex, double2& dir);
bool gjk_simplex(static_vector<double4, 4>& simplex, double4& dir);

// TODO use some kind of error measure instead of number of iterations
//
template<typename SupportA, typename SupportB>
inline int gjk_classify(const SupportA& support_a, const SupportB& support_b, double2& dir, uint max_iterations = 20) {
	static_vector<double2, 3> simplex;
	for (auto i : range(max_iterations)) {
		auto s = support_a(dir) - support_b(-dir);
		if (dot(s, dir) < 0)
			return 1;  // disjoint
		simplex.push_back(s);
		if (gjk_simplex(simplex, dir))
			return -1;  // overlapping
	}
	return 0;
}

template<typename SupportA, typename SupportB>
inline int gjk_classify(const SupportA& support_a, const SupportB& support_b, double4& dir, uint max_iterations = 20) {
	static_vector<double4, 4> simplex;
	for (auto i : range(max_iterations)) {
		auto s = support_a(dir) - support_b(-dir);
		if (dot(s, dir) < 0)
			return 1;  // disjoint
		simplex.push_back(s);
		if (gjk_simplex(simplex, dir))
			return -1;  // overlapping
	}
	return 0;
}

inline double4 sphere_support(double4 dir, double radius) {
	return dir * (radius / length(dir));
}

inline double4 capsule_support(double4 dir, double size, double radius) {
	if (dir.x > 0)
		return double4{size, 0, 0, 1} + dir * (radius / length(dir));
	return double4{-size, 0, 0, 1} + dir * (radius / length(dir));
}

inline double4 polyhedron_support(double4 dir, cspan<double4> vertices) {
	double md = -INF;
	double4 mv;
	for (double4 v : vertices) {
		double d = dot(v, dir);
		if (d > md) {
			md = d;
			mv = v;
		}
	}
	return mv;
}

inline double4 box_support(double4 dir, double4 size) {
	return sign_no_zero(dir) * size;
}
