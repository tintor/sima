#include <core/format.h>
#include <core/std.h>

#include <array>
#include <catch.hpp>

void compute(span<int> a) { a[0] = 2; }

TEST_CASE("span", "[span]") {
    vector<int> m = {0, 0, 0};
    m[1] = 0;
    span<int> s = {m.data(), m.size()};
    compute(s.last(2));
    REQUIRE(m[1] == 2);
}
