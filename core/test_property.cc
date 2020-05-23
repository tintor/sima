#include <core/property.h>

#include <catch.hpp>

struct M {
    long value;

    Property(M) {
        void operator=(const long i) { parent->value = i; }
        operator long() const { return parent->value; }
    } a;

    Property(M) {
        void operator=(const long i) { parent->value = i * 2; }
        operator long() const { return parent->value / 2; }
    } half_a;
};

static_assert(sizeof(M) == 8);

TEST_CASE("Property", "[property]") {
    M m;
    m.a = 10;
    REQUIRE(m.value == 10);
    REQUIRE(m.a == 10);
    m.half_a = 6;
    REQUIRE(m.value == 12);
    REQUIRE(m.half_a == 6);
}
