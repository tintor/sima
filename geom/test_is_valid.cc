#include <catch.hpp>
#include <geom/is_valid.h>
#include <core/string_util.h>

bool is_valid2(string_view text) {
	polygon2 poly;
	for (string_view t : split(text, ',')) {
		auto s = split(t, ' ');
		REQUIRE(s.size() == 2);
		poly.push_back(double2{parse<double>(s[0]), parse<double>(s[1])});
	}
	return IsValid(poly);
}

TEST_CASE("is_valid(polygon2)", "[is_valid]") {
	// not enough points
	REQUIRE(!is_valid2(""));
	REQUIRE(!is_valid2("0 0"));
	REQUIRE(!is_valid2("0 0, 0 1"));

	// valid triangle
	REQUIRE(is_valid2("0 0, 0 1, 1 0"));
	// valid triangle (opposite orientation)
	REQUIRE(is_valid2("0 0, 1 0, 0 1"));
	// valid square
	REQUIRE(is_valid2("0 0, 1 0, 1 1, 0 1"));

	// repeated point
	REQUIRE(!is_valid2("0 0, 0 0, 0 1, 1 0"));

	// self-overlaping edge
	REQUIRE(!is_valid2("0 0, 2 0, 1 0, 0 1"));

	// bow-tie without shared vertex
	REQUIRE(!is_valid2("0 0, 2 2, 2 0, 0 2"));

	// bow-tie with shared vertex
	REQUIRE(!is_valid2("0 0, 1 1, 2 2, 2 0, 1 1, 0 2"));

	// vertex touching edge inside
	REQUIRE(!is_valid2("0 0, 1 2, 2 0, 2 2, 0 2"));
	// edge touching edge inside
	REQUIRE(!is_valid2("0 0, 1 2, 1.1 2, 2 0, 2 2, 0 2"));
}
