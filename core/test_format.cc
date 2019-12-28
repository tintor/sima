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

	REQUIRE("00" == format("%02d", 0));
	REQUIRE("67" == format("%02d", 67));
	REQUIRE("000" == format("%03d", 0));
	REQUIRE("003" == format("%03d", 3));
	REQUIRE("034" == format("%03d", 34));

	REQUIRE("2" == format("%d", 2.1));
	REQUIRE("2.1" == format("%.1f", 2.1));

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

TEST_CASE("format - atomic", "[format]") {
	atomic<long> a;
	a = 3;
	REQUIRE("3" == format("%s", a));
}

TEST_CASE("format - int scientific", "[format]") {
	REQUIRE("0" == format("%h", 0));

	REQUIRE("9" == format("%h", 9));
	REQUIRE("99" == format("%h", 99));
	REQUIRE("1k" == format("%h", 999));

	REQUIRE("1k" == format("%h", 1000));
	REQUIRE("1m" == format("%h", 1000000));
	REQUIRE("1g" == format("%h", 1000000000lu));
	REQUIRE("1t" == format("%h", 1000000000000lu));

	REQUIRE("1400k" == format("%h", 1400000));
	REQUIRE("1235k" == format("%h", 1234567));
	REQUIRE("2g" == format("%h", 1999999999));
	REQUIRE("2m" == format("%h", 1999999));
	REQUIRE("2k" == format("%h", 1999));
	REQUIRE("2m" == format("%h", 1990009));

	REQUIRE("-9" == format("%h", -9));
	REQUIRE("-1k" == format("%h", -999));
	REQUIRE("-1400k" == format("%h", -1400000));
}

TEST_CASE("format - int time", "[format]") {
	REQUIRE("0" == format("%t", 0));
	REQUIRE("9" == format("%t", 9));
	REQUIRE("59" == format("%t", 59));
	REQUIRE("1:00" == format("%t", 60));
	REQUIRE("1:01" == format("%t", 61));
	REQUIRE("10:00" == format("%t", 600));
	REQUIRE("1:00:00" == format("%t", 3600));

	REQUIRE("-1:00:00" == format("%t", -3600));
	REQUIRE("-1:00" == format("%t", -60));
}

TEST_CASE("format - vector", "[format]") {
}
