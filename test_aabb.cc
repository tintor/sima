#include "catch.hpp"
#include "glm.h"
#include "aabb.h"
#include "format.h"
#include "util.h"
#include <random>

TEST_CASE("aabb basic", "[aabb]") {
	aabb<ivec3> box({ivec3{2, 7, 0}, ivec3{4, -1, 3}});
	REQUIRE(equal(box.min, ivec3{2, -1, 0}));
	REQUIRE(equal(box.max, ivec3{4, 7, 3}));
	REQUIRE(equal(box.size(), ivec3{2, 8, 3}));
	REQUIRE(equal(box.center(), ivec3{3, 3, 1}));
	REQUIRE(format("%s", box) == "aabb(2 -1 0, 4 7 3)");
}

TEST_CASE("aabb::operator==", "[aabb]") {
	ivec3 a = {0, 0, 0}, b = {1, 1, 1}, c = {1, 1, 0};
	aabb<ivec3> ab({a, b}), ac({a, c});
	REQUIRE(ab == ab);
	REQUIRE(!(ab != ab));
	REQUIRE(!(ab == ac));
	REQUIRE(ab != ac);
}

TEST_CASE("aabb::add()", "[aabb]") {
	ivec3 a = {0, 0, 0}, b = {1, 1, 1};
	aabb<ivec3> box({a, b});

	box.add(a);
	REQUIRE(equal(box.min, a));
	REQUIRE(equal(box.max, b));

	ivec3 c = {0, 0, 1}, d = {0, 4, 0};
	box.add(a - c);
	REQUIRE(equal(box.min, a - c));
	REQUIRE(equal(box.max, b));
	box.add(b + d);
	REQUIRE(equal(box.min, a - c));
	REQUIRE(equal(box.max, b + d));
}

template<typename RandomEngine>
ivec4 random_ivec4(RandomEngine& rnd, int a, int b) {
	std::uniform_int_distribution<int> dist(a, b);
	return ivec4{dist(rnd), dist(rnd), dist(rnd), dist(rnd)};
}

TEST_CASE("aabb(array_cptr<ivec4>) simd", "[aabb]") {
	std::default_random_engine rnd;
	std::vector<ivec4> w;
	w.resize(6);
	std::uniform_int_distribution<int> dist1(0, 1), dist2(1, 6);
	for (auto i : range(100)) {
		array_ptr<ivec4> v(w.data() + dist1(rnd), w.data() + dist2(rnd));
		if (v.size() == 0)
			continue;

		for (auto& e : v)
			e = random_ivec4(rnd, 0, 100);

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
}
