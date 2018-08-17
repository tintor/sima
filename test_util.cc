#include "util.h"
#include "string_util.h"
#include "catch.hpp"

TEST_CASE("min", "[util]") {
	REQUIRE(min(2, 1) == 1);
	REQUIRE(min(5, 7, 9) == 5);
	REQUIRE(min(9, 7, 5) == 5);
	REQUIRE(min(3, 0, -1, 6, 10, -4, 3) == -4);
}

TEST_CASE("split", "[util]") {
	using namespace std::literals;
	REQUIRE(split(""sv) == std::vector<std::string_view>{});
	REQUIRE(split("x"sv) == std::vector{"x"sv});
	REQUIRE(split(" a ana b[anana  "sv) == std::vector{"a"sv, "ana"sv, "b[anana"sv});

	REQUIRE(split(" an|a ba|na"sv) == std::vector{"an|a"sv, "ba|na"sv});
	REQUIRE(split(" an|a ba|na"sv, '|') == std::vector{" an"sv, "a ba"sv, "na"sv});
}
