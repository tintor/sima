#include <catch.hpp>
#include <core/exception.h>
#include <core/range.h>
#include <geom/pose.h>

TEST_CASE("pose basic", "[pose]") {
	pose2 a(d2(3, 1), PI / 5);
	pose2 b = inv(a);
	pose2 ab = mul(a, b);
	REQUIRE_NEAR(ab.position, d2(0, 0), 1e-15);
	REQUIRE_NEAR(ab.orientation, 0, 1e-15);
}

TEST_CASE("slerp angle", "[pose]") {
	double a = PI / 10;
	REQUIRE_NEAR(a / 4, slerp(0, a, 0.25), 1e-15);
	REQUIRE_NEAR(a - a / 4, slerp(a, 0, 0.25), 1e-15);

	REQUIRE_NEAR(PI - a / 2, slerp(PI - a, PI + a, 0.25), 1e-15);
	REQUIRE_NEAR(PI + a / 2, slerp(PI + a, PI - a, 0.25), 1e-15);
}
