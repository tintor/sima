#include <core/tensor.h>

#include <catch.hpp>

TEST_CASE("tensor_shape") {
    REQUIRE(tensor_shape().size() == 0);
    REQUIRE(tensor_shape(3).size() == 1);
    REQUIRE(tensor_shape(0, 1, 1).size() == 3);
    REQUIRE(tensor_shape(0, 1, 1, 4).size() == 4);

    REQUIRE(tensor_shape(1).volume() == 1);
    REQUIRE(tensor_shape(0, 2).volume() == 0);
    REQUIRE(tensor_shape(2, 3, 4, 5).volume() == 2*3*4*5);

    REQUIRE(tensor_shape(0, 3, 4, 5).remove_zeros() == tensor_shape(3, 4, 5));

    REQUIRE(tensor_shape(2) == tensor_shape(2));
    REQUIRE(tensor_shape(2, 2) == tensor_shape(2, 2));
    REQUIRE(tensor_shape(0, 2) != tensor_shape(2));
    REQUIRE(tensor_shape(2, 1) != tensor_shape(2));

    REQUIRE(tensor_shape(2, 1).pop_front() == tensor_shape(1));
    REQUIRE(tensor_shape(2).pop_front() == tensor_shape());
    REQUIRE(tensor_shape().pop_front() == tensor_shape());

    REQUIRE(tensor_shape(2, 1).pop_back() == tensor_shape(2));
    REQUIRE(tensor_shape(2).pop_back() == tensor_shape());
    REQUIRE(tensor_shape().pop_back() == tensor_shape());

    REQUIRE(tensor_shape(2, 1).push_front(3) == tensor_shape(3, 2, 1));
    REQUIRE(tensor_shape(2, 1).push_front(0) == tensor_shape(0, 2, 1));

    REQUIRE(tensor_shape(2, 4, 5).first(1) == tensor_shape(2));
    REQUIRE(tensor_shape(2, 4, 5).first(2) == tensor_shape(2, 4));
    REQUIRE(tensor_shape(0, 4, 5).first(3) == tensor_shape(0, 4, 5));

    REQUIRE(tensor_shape(2, 4, 5).last(1) == tensor_shape(5));
    REQUIRE(tensor_shape(2, 4, 5).last(2) == tensor_shape(4, 5));
    REQUIRE(tensor_shape(0, 4, 5).last(3) == tensor_shape(0, 4, 5));

    REQUIRE(string(tensor_shape(2, 4, 5)) == "[2 4 5]");
}

TEST_CASE("tensor") {
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
    REQUIRE(tensor(m.data(), {3}) != tensor(n.data(), {3, 1}));

    REQUIRE(tensor(m.data(), {2, 3}).slice(0) == tensor(m.data(), {3}));
    REQUIRE(tensor(m.data(), {2, 3}).slice(1) == tensor(m.data() + 3, {3}));

    REQUIRE(tensor(m.data(), {2}));
    REQUIRE(!tensor());
    // REQUIRE(!tensor(nullptr, {}));

    REQUIRE(tensor().size() == 0);
    REQUIRE(tensor() == tensor());
    REQUIRE(tensor(m.data(), {1}) != tensor());
}

TEST_CASE("vtensor") {
    vector<int> a;
    REQUIRE(a.data() == nullptr);
}
