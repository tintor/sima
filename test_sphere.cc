#include "sphere.h"
#include "range.h"
#include "catch.hpp"
#include "aabb.h"

TEST_CASE("minimal_sphere(sphere, sphere)", "[sphere]") {
	std::default_random_engine rnd;
	double e = 1;
	while (e > 1e-10) {
		sphere a = make_sphere(uniform3(rnd, -1, 1) * 1000.0, 1);
		sphere b = make_sphere(center(a) + e * uniform_dir3(rnd), 1);
		sphere c = minimal_sphere(a, b);
		REQUIRE(radius(c) >= radius(a));
		REQUIRE(radius(c) >= radius(b));
		REQUIRE(squared(center(c) - center(a)) <= squared(radius(c) - radius(a)) * 1.01);
		REQUIRE(squared(center(c) - center(b)) <= squared(radius(c) - radius(b)) * 1.01);
		e *= 0.9;
	}
}

sphere minimal_sphere_brute_force(span<const point3> points) {
	if (points.size() == 1)
		return make_sphere(points[0], 0);
	if (points.size() == 2)
		return minimal_sphere(points[0], points[1]);
	if (points.size() == 3)
		return minimal_sphere(points[0], points[1], points[2]);

	sphere minimal = make_sphere({0, 0, 0}, -1);
	for (auto a : range(points.size()))
		for (auto b : range(a + 1, points.size()))
			for (auto c : range(b + 1, points.size()))
				for (auto d : range(c + 1, points.size())) {
					sphere s = minimal_sphere(points[a], points[b], points[c], points[d]);
					if (radius(s) > radius(minimal))
						minimal = s;
				}
	return minimal;
}

double max_distance(point3 center, span<const point3> points) {
	double d2 = 0;
	for (point3 p : points)
		d2 = max(d2, squared(p - center));
	return sqrt(d2);
}

bool contains(sphere s, span<const point3> points) {
	for (auto p : points)
		if (!contains(s, p))
			return false;
	return true;
}

sphere bounding_sphere_iterative(span<const point3> points, double eps) {
	// first approximation is center of aabb
	aabb3p box(points);
	point3 center = box.center();
	double radius = max_distance(center, points);

	double max_shift = radius * 0.5;
	std::default_random_engine rnd;
	while (max_shift > eps) {
		point3 new_center = center + uniform3(rnd, -max_shift, max_shift);
		double new_radius = max_distance(new_center, points);
		if (new_radius < radius) {
			center = new_center;
			radius = new_radius;
			max_shift = min(max_shift, radius * 0.5);
		} else {
			max_shift *= 0.99;
		}
	}
	sphere s = make_sphere(center, radius);
	while (!contains(s, points)) {
		radius = std::nexttoward(radius, std::numeric_limits<double>::max());
		s = make_sphere(center, radius);
	}
	return s;
}

TEST_CASE("minimal_sphere(point3 x2)", "[sphere]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	for (auto i : range(20000)) {
		point3 a = uniform3(rnd, -1, 1);
		point3 b = uniform3(rnd, -1, 1);
		sphere s1 = minimal_sphere(a, b);
		sphere s2 = bounding_sphere_iterative({a, b}, 1e-10);

		REQUIRE(contains(s1, a));
		REQUIRE(contains(s1, b));

		REQUIRE(contains(s2, a));
		REQUIRE(contains(s2, b));

		REQUIRE(radius(s1) <= radius(s2) + 1e-15);
	}
}

TEST_CASE("minimal_sphere(point3 x3)", "[sphere]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	for (auto i : range(20000)) {
		point3 a = uniform3(rnd, -1, 1);
		point3 b = uniform3(rnd, -1, 1);
		point3 c = uniform3(rnd, -1, 1);
		sphere s1 = minimal_sphere(a, b, c);
		sphere s2 = bounding_sphere_iterative({a, b, c}, 1e-10);

		REQUIRE(contains(s1, a));
		REQUIRE(contains(s1, b));
		REQUIRE(contains(s1, c));

		REQUIRE(contains(s2, a));
		REQUIRE(contains(s2, b));
		REQUIRE(contains(s2, c));

		REQUIRE(radius(s1) <= radius(s2) + 1e-15);
	}
}

TEST_CASE("minimal_sphere(point3 x4)", "[sphere]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	for (auto i : range(20000)) {
		point3 a = uniform3(rnd, -1, 1);
		point3 b = uniform3(rnd, -1, 1);
		point3 c = uniform3(rnd, -1, 1);
		point3 d = uniform3(rnd, -1, 1);
		sphere s1 = minimal_sphere(a, b, c, d);
		sphere s2 = bounding_sphere_iterative({a, b, c, d}, 1e-10);

		REQUIRE(contains(s1, a));
		REQUIRE(contains(s1, b));
		REQUIRE(contains(s1, c));
		REQUIRE(contains(s1, d));

		REQUIRE(contains(s2, a));
		REQUIRE(contains(s2, b));
		REQUIRE(contains(s2, c));
		REQUIRE(contains(s2, d));

		REQUIRE(radius(s1) <= radius(s2) + 1e-15);
	}
}

TEST_CASE("minimal_sphere(span<point3>)", "[sphere]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	aligned_vector<point3> points;
	for (auto i : range(5000)) {
		std::uniform_int_distribution<int> dist(5, 16);
		int size = dist(rnd);

		points.resize(size);
		for (auto j : range(size))
			points[j] = uniform3(rnd, -1.0, 1.0);

		sphere s2 = minimal_sphere_brute_force(points);
		sphere s1 = minimal_sphere(points);
		REQUIRE(abs(radius(s1) - radius(s2)) <= 1e-10);
	}
}

TEST_CASE("bounding_sphere(span<point3>)_bound", "[sphere][!hide]]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	aligned_vector<point3> points;
	double4 e;
	for (auto i : range(50000)) {
		int size = 1000;

		points.resize(size);
		for (auto j : range(size))
			points[j] = uniform3(rnd, -1.0, 1.0);

		e += bounding_sphere(points);
	}
	print("e=%s\n", e);
}

TEST_CASE("bounding_sphere(span<point3>)_min", "[sphere][!hide]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	aligned_vector<point3> points;
	double4 e;
	for (auto i : range(50000)) {
		int size = 1000;

		points.resize(size);
		for (auto j : range(size))
			points[j] = uniform3(rnd, -1.0, 1.0);

		e += minimal_sphere(points);
	}
	print("e=%s\n", e);
}

TEST_CASE("bounding_sphere(span<point3>)", "[sphere]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	aligned_vector<point3> points;
	for (auto i : range(25000)) {
		int size = 100;

		points.resize(size);
		for (auto j : range(size))
			points[j] = uniform3(rnd, -1.0, 1.0);

		sphere s1 = minimal_sphere(points);
		sphere s2 = bounding_sphere(points);
		REQUIRE(abs(radius(s1) - radius(s2)) <= 0.05);
	}
}
