#include <core/string_util.h>
#include <catch.hpp>

TEST_CASE("split", "[util]") {
	using namespace std::literals;
	REQUIRE(split(""sv) == vector<string_view>{});
	REQUIRE(split("x"sv) == vector{"x"sv});
	REQUIRE(split(" a ana b[anana  "sv) == vector{"a"sv, "ana"sv, "b[anana"sv});

	REQUIRE(split(" an|a ba|na"sv) == vector{"an|a"sv, "ba|na"sv});
	REQUIRE(split(" an|a ba|na"sv, '|') == vector{" an"sv, "a ba"sv, "na"sv});
}

TEST_CASE("natural_less", "[string_util]") {
	using namespace std::literals;
	REQUIRE(natural_less("ma2", "ma10"));
	REQUIRE(!natural_less("ma10", "ma2"));
	REQUIRE(!natural_less("ma2", "ma2"));
	REQUIRE(natural_less("a", "b"));
}

TEST_CASE("cat", "[string_util]") {
	using namespace std::literals;
	REQUIRE("abcde" == cat("abc"sv, "de"s));
	REQUIRE("abcde" == cat("abc"sv, "de"sv));
}
