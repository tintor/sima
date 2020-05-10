#include <core/array_deque.h>
#include <core/range.h>

#include <catch.hpp>

TEST_CASE("array_deque basic", "[array_deque]") {
    array_deque<int> deque;
    REQUIRE(deque.capacity() == 0);

    for (int i : range(100)) {
        REQUIRE(deque.size() == i);
        deque.push_back(i);
    }
    REQUIRE(deque.capacity() == 128);

    for (int i : range(40)) REQUIRE(i == deque.pop_front());

    REQUIRE(deque.size() == 60);
    REQUIRE(deque.capacity() == 128);

    for (int i : range(100, 170)) deque.push_back(i);
    REQUIRE(deque.size() == 130);

    for (int i : range(130)) REQUIRE(deque[i] == 40 + i);
}
