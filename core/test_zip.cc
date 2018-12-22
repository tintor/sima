#include <core/zip.h>
#include <core/std.h>
#include "catch.hpp"

void test_zip(cspan<int> a, cspan<int> b, cspan<pair<int, int>> exp) {
	vector<pair<int, int>> c;
	for (auto e : czip(a, b))
		c.push_back(e);
	REQUIRE(cspan<pair<int, int>>(c) == exp);
}

TEST_CASE("zip") {
	test_zip({1, 2, 3}, {}, {});
	test_zip({}, {1, 2, 3}, {});
	test_zip({1, 2, 3}, {4, 5}, {{1, 4}, {2, 5}});
	test_zip({1, 2}, {3, 4, 5}, {{1, 3}, {2, 4}});
}
