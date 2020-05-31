#include <geom/vector.h>

#include <catch.hpp>
#include <klein/klein.hpp>

using namespace kln;

TEST_CASE("geom_algebra basic", "[geom_algebra]") {
    plane p1{1.f, 2.f, 3.f, 4.f};
    plane p2{2.f, 3.f, -1.f, -2.f};
    motor p12 = p1 * p2;
}

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
