#include "catch.hpp"
#include "glm.h"
#include "aabb.h"
#include "format.h"

TEST_CASE("aabb basic", "[aabb]") {
	aabb<ivec3> box({ivec3(2, 7, 0), ivec3(4, -1, 3)});
	REQUIRE(box.min == ivec3(2, -1, 0));
	REQUIRE(box.max == ivec3(4, 7, 3));
	REQUIRE(box.size() == ivec3(2, 8, 3));
	REQUIRE(box.center() == ivec3(3, 3, 1));
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
	REQUIRE(box.min == a);
	REQUIRE(box.max == b);

	ivec3 c = {0, 0, 1}, d = {0, 4, 0};
	box.add(a - c);
	REQUIRE(box.min == a - c);
	REQUIRE(box.max == b);
	box.add(b + d);
	REQUIRE(box.min == a - c);
	REQUIRE(box.max == b + d);
}
