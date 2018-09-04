#include "generators.h"
#include "properties.h"
#include "format.h"
#include "catch.hpp"
#include <random>

TEST_CASE("generate_cylinder", "[generators]") {
	double r = 100;
	double h = 20;
	auto m = generate_cylinder(500, r, r, -h/2, h/2);

	double ev = r * r * h * M_PI;
	double av = Volume(m);
	REQUIRE(abs(av - ev) / ev < 0.01);
}

TEST_CASE("generate_cone", "[generators]") {
	double r = 100;
	double h = 20;
	auto m = generate_cylinder(500, r, 0, -h/2, h/2);

	double ev = r * r * h * M_PI / 3;
	double av = Volume(m);
	REQUIRE(abs(av - ev) / ev < 0.01);
}

TEST_CASE("generate_sphere", "[generators]") {
	std::default_random_engine rnd;
	int n = 500;
	double r = 1000.0;
	auto m = generate_sphere(n, r, rnd);

	double ev = 4.0 / 3 * M_PI * r * r * r;
	double av = Volume(m);
	REQUIRE(m.size() == n * 2 - 4);
	REQUIRE(abs(av - ev) / ev < 0.012);
}

TEST_CASE("generate_regular_polyhedra", "[generators][.]") {
	// tetrahedron
	//generate_regular_polyhedra({{0, 1, 2}, {3, 1, 0}, {3, 2, 1}, {3, 0, 2}});
	// tetrahedron
	generate_regular_polyhedra2({{3, 4}});
}
