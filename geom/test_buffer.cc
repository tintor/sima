#include <catch.hpp>
#include <geom/buffer.h>
#include <core/range.h>
#include <geom/polygon.h>
#include <geom/properties.h>

static double length(cspan<double2> chain) {
	double m = 0;
	const int n = chain.size();
	for (int i : range(n))
		m += length(chain[i] - chain[(i + 1) % n]);
	return m;
}

static double shortest_edge(cspan<double2> polygon) {
	double m = INF;
	const int n = polygon.size();
	for (int i : range(n))
		minimize(m, squared(polygon[i] - polygon[(i + 1) % n]));
	return sqrt(m);
}

constexpr int Vertices = 128;

TEST_CASE("buffer point", "[buffer]") {
	array<double2, 1> s = { double2{2, 4} };
	const double r = 2;
	auto b = ComputeBuffer(s, r, Vertices);
	REQUIRE(area(b) == Approx(r * r * PI).epsilon(1e-3));
	REQUIRE(length(b) == Approx(2 * r * PI).epsilon(1e-3));
	REQUIRE(length(centroid(b) - s[0]) <= 1e-6);
	REQUIRE(shortest_edge(b) > 2 * r * PI / (Vertices + 1));
}

TEST_CASE("buffer line", "[buffer]") {
	array<double2, 2> s = { double2{2, 0}, double2{3, 0} };
	const double r = 2;
	auto b = ComputeBuffer(s, r, Vertices);
	double ea = r * r * PI + length(s) * r;
	REQUIRE(area(b) == Approx(ea).epsilon(1e-3));
	REQUIRE(length(b) == Approx(length(s) + 2 * r * PI).epsilon(1e-3));
	REQUIRE(length(centroid(b) - (s[0] + s[1]) / 2) <= 1e-6);
	REQUIRE(shortest_edge(b) > 2 * r * PI / (Vertices + 1));
}

TEST_CASE("buffer triangle", "[buffer]") {
	array<double2, 3> s = { double2{0, 0}, double2{2, 0}, double2{0, 1} };
	const double r = 1;
	auto b = ComputeBuffer(s, r, Vertices);
	double ea = area(s) + r * r * PI + length(s) * r;
	REQUIRE(area(b) == Approx(ea).epsilon(1e-3));
	REQUIRE(length(b) == Approx(length(s) + 2 * r * PI).epsilon(1e-3));
	REQUIRE(shortest_edge(b) > 2 * r * PI / (Vertices * 1.2));
}
