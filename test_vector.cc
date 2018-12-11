#include "vector.h"
#include "catch.hpp"

TEST_CASE("all", "[vector]") {
	REQUIRE(all(int4{-1, -1, -1, -1}));
	REQUIRE(!all(int4{0, -1, -1, -1}));
	REQUIRE(!all(int4{-1, 0, -1, -1}));
	REQUIRE(!all(int4{-1, -1, 0, -1}));
	REQUIRE(!all(int4{-1, -1, -1, 0}));
}

TEST_CASE("any", "[vector]") {
	REQUIRE(!any(int4{0, 0, 0, 0}));
	REQUIRE(any(int4{-1, 0, 0, 0}));
	REQUIRE(any(int4{0, -1, 0, 0}));
	REQUIRE(any(int4{0, 0, -1, 0}));
	REQUIRE(any(int4{0, 0, 0, -1}));
	REQUIRE(any(int4{-1, -1, -1, -1}));
}

TEST_CASE("sign", "[vector]") {
	double4 a = sign_no_zero(double4{0.0, 2.1, -0.14, -1e100});
	REQUIRE(a.x == 1.0);
	REQUIRE(a.y == 1.0);
	REQUIRE(a.z == -1.0);
	REQUIRE(a.w == -1.0);
	REQUIRE(all(a == double4{1.0, 1.0, -1.0, -1.0}));
}
