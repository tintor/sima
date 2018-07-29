#include "auto.h"
#include "catch.hpp"

TEST_CASE("auto") {
	int a = 1;
	{
		ON_SCOPE_EXIT(a *= 2);
		a += 10;
	}
	REQUIRE(a == 22);
}
