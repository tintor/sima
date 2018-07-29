#include "array_set.h"
#include "range.h"
#include "catch.hpp"

TEST_CASE("array_set") {
	array_set<int> a;
	REQUIRE(a.insert(3));
	REQUIRE(a.insert(7));
	REQUIRE(!a.insert(3));
	REQUIRE(a.size() == 2);
	REQUIRE(a.contains(3));
	REQUIRE(a.contains(7));
	for (auto i : range(1, 100))
		if (i != 3 && i != 7)
			REQUIRE(!a.contains(i));
}

TEST_CASE("array_set big") {
	array_set<int> a;
	for (int e : range(1, 1001)) {
		REQUIRE(a.insert(e));
		a.print();
	}
	REQUIRE(a.size() == 1000);
	for (int e : range(1, 1001))
		REQUIRE(a.erase(e));
	REQUIRE(a.size() == 0);
}

TEST_CASE("array_map big") {
	array_map<int, float> a;
	for (int e : range(1, 1001)) {
		REQUIRE(a.insert(e, e * 0.1f));
		a.print();
	}
	REQUIRE(a.size() == 1000);
	for (int e : range(1, 1001))
		REQUIRE(a.erase(e));
	REQUIRE(a.size() == 0);
}
