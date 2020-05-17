#include <core/diff.h>

#include <catch.hpp>

TEST_CASE("diff") {
    auto c = Constant(5);
    auto relu = Relu(c);
    auto model = TopoSort({relu});
    REQUIRE(model == vector<PDiff>{c, relu});
    Setup(model, 1);
    Forward(model);
    REQUIRE(relu->v.shape() == tensor_shape{1});
    REQUIRE(relu->v[0] == 5);
    Backward(model);

//    PDiff data = Data({0, 2});
//    ZeroInit zinit;
//    PDiff param = Parameter({2}, zinit);
}
