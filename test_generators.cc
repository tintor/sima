#include "generators.h"
#include "properties.h"
#include "format.h"
#include "catch.hpp"
#include <random>

TEST_CASE("generate_cylinder") {
	double r = 100;
	double h = 20;
	auto m = generate_cylinder(500, r, r, -h/2, h/2);

	double ev = r * r * h * M_PI;
	double av = volume(m);
	REQUIRE(abs(av - ev) / ev < 0.01);
}

TEST_CASE("generate_cone") {
	double r = 100;
	double h = 20;
	auto m = generate_cylinder(500, r, 0, -h/2, h/2);

	double ev = r * r * h * M_PI / 3;
	double av = volume(m);
	REQUIRE(abs(av - ev) / ev < 0.01);
}

TEST_CASE("generate_sphere") {
	std::default_random_engine rnd;
	int n = 500;
	double r = 1000.0;
	auto m = generate_sphere(n, r, rnd);

	double ev = 4.0 / 3 * M_PI * r * r * r;
	double av = volume(m);
	REQUIRE(m.size() == n * 2 - 4);
	REQUIRE(abs(av - ev) / ev < 0.012);
}
