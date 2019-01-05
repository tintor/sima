#include <core/std.h>
#include <core/range.h>
#include <core/format.h>
#include <core/dynamic_array.h>
#include <geom/vector.h>

bool gjk_simplex(static_vector<double2, 3>& simplex, double2& dir);
bool gjk_simplex(static_vector<double3, 4>& simplex, double3& dir);

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
inline int gjk_classify(const SupportA& support_a, const SupportB& support_b, double3& dir, uint max_iterations = 20) {
	static_vector<double3, 4> simplex;
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

inline double3 sphere_support(double3 dir, double radius) {
	return dir * (radius / length(dir));
}

inline double3 capsule_support(double3 dir, double size, double radius) {
	if (dir.x > 0)
		return double3{size, 0, 0} + dir * (radius / length(dir));
	return double3{-size, 0, 0} + dir * (radius / length(dir));
}

inline double3 polyhedron_support(double3 dir, cspan<double3> vertices) {
	double md = -INF;
	double3 mv;
	for (double3 v : vertices) {
		double d = dot(v, dir);
		if (d > md) {
			md = d;
			mv = v;
		}
	}
	return mv;
}

inline double3 box_support(double3 dir, double3 size) {
	return sign_no_zero(dir) * size;
}
