#include "sphere.h"
#include <core/range.h>
#include "catch.hpp"
#include <geom/aabb.h>

TEST_CASE("minimal_sphere(sphere, sphere)", "[sphere]") {
	std::default_random_engine rnd;
	double e = 1;
	while (e > 1e-10) {
		sphere a = sphere(uniform3(rnd, -1, 1) * 1000.0, 1);
		sphere b = sphere(a.center() + e * uniform_dir3(rnd), 1);
		sphere c = minimal_sphere(a, b);
		REQUIRE(c.radius() >= a.radius());
		REQUIRE(c.radius() >= b.radius());
		REQUIRE(squared(c.center() - a.center()) <= squared(c.radius() - a.radius()) * 1.01);
		REQUIRE(squared(c.center() - b.center()) <= squared(c.radius() - b.radius()) * 1.01);
		e *= 0.9;
	}
}

sphere minimal_sphere_brute_force(cspan<point3> points) {
	if (points.size() == 1)
		return sphere(points[0], 0);
	if (points.size() == 2)
		return minimal_sphere(points[0], points[1]);
	if (points.size() == 3)
		return minimal_sphere(points[0], points[1], points[2]);

	sphere minimal = sphere(double4{0, 0, 0, 1}, -1);
	for (auto a : range(points.size()))
		for (auto b : range(a + 1, points.size()))
			for (auto c : range(b + 1, points.size()))
				for (auto d : range(c + 1, points.size())) {
					sphere s = minimal_sphere(points[a], points[b], points[c], points[d]);
					if (s.radius() > minimal.radius())
						minimal = s;
				}
	return minimal;
}

double max_distance(point3 center, cspan<point3> points) {
	double d2 = 0;
	for (point3 p : points)
		d2 = max(d2, squared(p - center));
	return sqrt(d2);
}

bool contains(sphere s, cspan<point3> points) {
	for (auto p : points)
		if (!s.contains(p))
			return false;
	return true;
}

sphere bounding_sphere_iterative(cspan<point3> points, double eps) {
	// first approximation is center of aabb
	aabb4 box(points);
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
	sphere s = sphere(center, radius);
	while (!contains(s, points)) {
		radius = std::nexttoward(radius, std::numeric_limits<double>::max());
		s = sphere(center, radius);
	}
	return s;
}

bool contains(sphere s, double4 a) { return s.contains(a); }
double radius(sphere s) { return s.radius(); }

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

		e += bounding_sphere(points).center();
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

		e += minimal_sphere(points).center();
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
