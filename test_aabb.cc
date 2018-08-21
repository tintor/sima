#include "catch.hpp"
#include "aabb.h"
#include "format.h"
#include "util.h"
#include <random>

TEST_CASE("aabb3 basic", "[aabb]") {
	aabb3 box(double3{2, 7, 0}, double3{4, -1, 3});
	REQUIRE(equal(box.min, double3{2, -1, 0}));
	REQUIRE(equal(box.max, double3{4, 7, 3}));
	REQUIRE(equal(box.size(), double3{2, 8, 3}));
	REQUIRE(equal(box.center(), double3{3, 3, 1.5}));
	REQUIRE(format("%.1lf", box) == "aabb(2.0 -1.0 0.0, 4.0 7.0 3.0)");
	REQUIRE(format("%.0lf", box) == "aabb(2 -1 0, 4 7 3)");
}

TEST_CASE("aabb::operator==", "[aabb]") {
	double3 a = {0, 0, 0}, b = {1, 1, 1}, c = {1, 1, 0};
	aabb3 ab(segment3{a, b}), ac(segment3{a, c});
	REQUIRE(ab == ab);
	REQUIRE(!(ab != ab));
	REQUIRE(!(ab == ac));
	REQUIRE(ab != ac);
}

TEST_CASE("aabb::add()", "[aabb]") {
	double3 a = {0, 0, 0}, b = {1, 1, 1};
	aabb3 box(segment3{a, b});

	box.add(a);
	REQUIRE(equal(box.min, a));
	REQUIRE(equal(box.max, b));

	double3 c = {0, 0, 1}, d = {0, 4, 0};
	box.add(a - c);
	REQUIRE(equal(box.min, a - c));
	REQUIRE(equal(box.max, b));
	box.add(b + d);
	REQUIRE(equal(box.min, a - c));
	REQUIRE(equal(box.max, b + d));
}

/*TEST_CASE("aabb(array_cptr<ivec4>) simd", "[aabb]") {
	std::default_random_engine rnd;
	std::vector<ivec4> w;
	w.resize(6);
	std::uniform_int_distribution<int> dist1(0, 1), dist2(1, 6);
	for (auto i : range(100)) {
		array_ptr<ivec4> v(w.data() + dist1(rnd), w.data() + dist2(rnd));
		if (v.size() == 0)
			continue;

		for (auto& e : v)
			e = uniform4(rnd, 0, 100);

		aabb<ivec4> box(v);
		ivec4 emin = v[0], emax = v[0];
		for (auto e : v) {
			emin.x = min(emin.x, e.x);
			emin.y = min(emin.y, e.y);
			emin.z = min(emin.z, e.z);
			emin.w = min(emin.w, e.w);

			emax.x = max(emax.x, e.x);
			emax.y = max(emax.y, e.y);
			emax.z = max(emax.z, e.z);
			emax.w = max(emax.w, e.w);
		}

		REQUIRE(equal(box.min, emin));
		REQUIRE(equal(box.max, emax));
	}
}*/
