#include <core/dynamic_array.h>
#include <core/std.h>
#include "catch.hpp"

TEST_CASE("dynamic_array") {
	dynamic_array<int> m;
	REQUIRE(m.begin() == nullptr);
	REQUIRE(m.size() == 0);

	m.resize(3);
	REQUIRE(m.begin() != nullptr);
	REQUIRE(m.size() == 3);
	m[0] = 10;
	m[1] = 20;
	m[2] = 30;
	vector<int> v;
	for (int e : m)
		v.push_back(e);
	REQUIRE(v == vector{10, 20, 30});

	m.resize(4);
	REQUIRE(m.size() == 4);
	m[3] = 30;
	REQUIRE(m[0] == 10);

	dynamic_array<int> n = std::move(m);
	REQUIRE(n.size() == 4);
	REQUIRE(m.size() == 0);
}
