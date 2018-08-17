#include "span.h"
#include "catch.hpp"
#include <array>

#include "format.h"
#include "triangle.h"

void compute(span<int> a) {
	a[0] = 2;
}

TEST_CASE("span") {
	std::vector<int> m;
	m[1] = 0;
	span<int> s = {m.data(), m.size()};
	compute(s.last(2));
	REQUIRE(m[1] == 2);
}
