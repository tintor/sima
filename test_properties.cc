#include "properties.h"
#include "generators.h"
#include "catch.hpp"

// TODO center of mass
// TODO translate mesh
// TODO rotate mesh
TEST_CASE("volume_of_cube") {
	auto m = generate_box(2, 3, 5);
	REQUIRE(signed_volume(m) == 240);
}
