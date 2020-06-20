#include <core/format.h>
#include <core/tensor.h>

#include <catch.hpp>

TEST_CASE("dim4", "[tensor]") {
    REQUIRE(dim4{}.rank() == 0);
    REQUIRE(dim4{1}.rank() == 0);
    REQUIRE(dim4{3}.rank() == 1);
    REQUIRE(dim4{1, 3, 1}.rank() == 1);
    REQUIRE(dim4{0, 3, 1}.rank() == -1);

    REQUIRE(dim4{1}.elements() == 1);
    REQUIRE(dim4{2, 3, 4, 5}.elements() == 2 * 3 * 4 * 5);

    REQUIRE(dim4{2} == dim4{2});
    REQUIRE(dim4{2, 2} == dim4{2, 2});
    REQUIRE(dim4{2, 1} == dim4{2});

    REQUIRE(dim4{2, 1}.pop_front() == dim4{1});
    REQUIRE(dim4{2}.pop_front() == dim4{});
    REQUIRE(dim4{}.pop_front() == dim4{});

    REQUIRE(dim4{2, 1}.pop_back() != dim4{2});
    REQUIRE(dim4{2}.pop_back() == dim4{});
    REQUIRE(dim4{}.pop_back() == dim4{});

    REQUIRE(dim4{2, 1}.push_front(3) == dim4{3, 2, 1});

    // REQUIRE(dim4{2, 4, 5}.first(1) == dim4{2});
    // REQUIRE(dim4{2, 4, 5}.first(2) == dim4{2, 4});

    // REQUIRE(dim4{2, 4, 5}.last(1) == dim4{5});
    // REQUIRE(dim4{2, 4, 5}.last(2) == dim4{4, 5});

    REQUIRE(dim4{2, 4, 5}.str() == "[2 4 5]");
    REQUIRE(dim4{2}.str() == "[2]");
    REQUIRE(dim4{1}.str() == "scalar");
    REQUIRE(dim4{}.str() == "scalar");
}

TEST_CASE("tensor", "[tensor]") {
    array<float, 6> m = {1, 2, 3, 4, 5, 6};
    array<float, 3> n = {1, 2, 3};
    REQUIRE(tensor(m.data(), {6})[2] == 3);

    REQUIRE(tensor(m.data(), {6})(3) == 4);
    REQUIRE(tensor(m.data(), {2, 3})(0, 0) == 1);
    REQUIRE(tensor(m.data(), {2, 3})(0, 1) == 2);
    REQUIRE(tensor(m.data(), {2, 3})(0, 2) == 3);
    REQUIRE(tensor(m.data(), {2, 3})(1, 0) == 4);
    REQUIRE(tensor(m.data(), {2, 3})(1, 1) == 5);
    REQUIRE(tensor(m.data(), {2, 3})(1, 2) == 6);

    REQUIRE(tensor(m.data(), {3}) == tensor(n.data(), {3}));
    REQUIRE(tensor(m.data() + 1, {3}) != tensor(n.data(), {3}));
    REQUIRE(tensor(m.data(), {3}) != tensor(n.data(), {2}));
    REQUIRE(tensor(m.data(), {3}) == tensor(n.data(), {3, 1}));

    REQUIRE(tensor(m.data(), {2, 3}).slice(0) == tensor(m.data(), {3}));
    REQUIRE(tensor(m.data(), {2, 3}).slice(1) == tensor(m.data() + 3, {3}));

    REQUIRE(tensor(m.data(), {2}));
    REQUIRE(!tensor());
    REQUIRE(!tensor(nullptr, {0}));

    REQUIRE(tensor().size == 0);
    REQUIRE(tensor() == tensor());
    REQUIRE(tensor(m.data(), {1}) != tensor());
}

TEST_CASE("vtensor", "[tensor]") {
    vector<int> a;
    REQUIRE(a.data() == nullptr);
}
