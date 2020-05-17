#include <santorini/model.h>

#include <catch.hpp>

TEST_CASE("santorini model") {
    ValueFunction model(100);
    model.Compile();
    model.Print();
}
