#include "zip.h"
#include "catch.hpp"

void test_zip(const std::vector<int>& a, const std::vector<int>& b, const std::vector<std::pair<int, int>>& exp) {
	std::vector<std::pair<int, int>> c;
	for (auto e : czip(a, b))
		c.push_back(e);
	REQUIRE(c == exp);
}

TEST_CASE("zip") {
	test_zip({1, 2, 3}, {}, {});
	test_zip({}, {1, 2, 3}, {});
	test_zip({1, 2, 3}, {4, 5}, {{1, 4}, {2, 5}});
	test_zip({1, 2}, {3, 4, 5}, {{1, 3}, {2, 4}});
}
