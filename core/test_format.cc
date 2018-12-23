#include <core/format.h>
#include <core/int.h>
#include <catch.hpp>
#include <iostream>

using namespace std::literals;

TEST_CASE("format - basic", "[format]") {
	REQUIRE("0" == format("%s", __int128(0)));
	REQUIRE("-1234567890" == format("%s", __int128(-1234567890)));

	string m;
	m << "Hello" << 10 << "world!"sv;
	REQUIRE("Hello10world!" == m);

	REQUIRE("0" == format("%s", int32_t(0)));
	REQUIRE("-1" == format("%s", int32_t(-1)));
	REQUIRE("-123456789" == format("%s", int32_t(-123456789)));
	REQUIRE("1" == format("%s", int32_t(1)));
	REQUIRE("123456789" == format("%s", int32_t(123456789)));

	REQUIRE("3.00" == format("%.2f", 3.0f));
	REQUIRE("3.000000" == format("%f", 3.0f));
	REQUIRE("3" == format("%s", 3.0f));
	REQUIRE("3.1" == format("%s", 3.1f));
	REQUIRE("3.14" == format("%s", 3.14f));

	REQUIRE("" == format(""));
	REQUIRE_THROWS_WITH(format("%s"), "format: not enough arguments");
	REQUIRE("%" == format("%%"));
	REQUIRE("e" == format("e"));
	REQUIRE("e%" == format("e%%"));
	REQUIRE("e1b" == format("e%sb", 1));
	REQUIRE("e1b2.3ctruea" == format("e%sb%sc%sa", 1, 2.3, true));
}

TEST_CASE("format - vector", "[format]") {
}
