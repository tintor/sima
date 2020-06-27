#include <core/std.h>
#include <core/thread.h>

#include <catch.hpp>

TEST_CASE("parallel_map_reduce", "[thread]") {
    std::atomic<int> next = 0;
    int result = parallel_map_reduce<int>([&next](){ int a = ++next; return a; }, [](int& acc, int e) { acc += e; });
    int c = next.load();
    REQUIRE(c * (c + 1) / 2 == result);
}
