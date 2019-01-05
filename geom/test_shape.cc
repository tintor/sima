#include <geom/shape.h>
#include <catch.hpp>

TEST_CASE("is_valid") {
/*	REQUIRE(is_valid(vector<triangle3>{}) == Validity::TooFewFaces);
	double3 a(0, 0, 0);
	double3 b(1, 0, 0);
	double3 c(0, 1, 0);
	double3 d(0, 0, 1);*/
}

TEST_CASE("SphereMesh - ConvexHull - IsConvex") {
/*	std::default_random_engine rnd;
	Mesh3d mesh = SphereMesh(100, rnd);
	real v = 4.0 / 3 * PI;
	REQUIRE(abs(volume(mesh) - v) / v < 0.15);
	REQUIRE(is_valid(mesh) == Validity::OK);
	REQUIRE(is_convex(mesh));*/
}

TEST_CASE("CreateCrossMesh") {
	//Mesh3d cross = CreateCrossMesh(0.05, 0.25);
	//REQUIRE(IsValid(cross) == Validity::OK);
}
