#include "tesselate.h"
#include "format.h"
#include "catch.hpp"
#include <unordered_set>
#include <random>
#include <iostream>

TEST_CASE("tesselate 2") {
	ivec2 a = {0, 0};
	ivec2 b = {2, 1};
	ivec2 c = {0, 2};
	ivec2 d = {1, 1};
	REQUIRE(tesselate({a, b, c, d}) == imesh2{{d, a, b}, {d, b, c}});
}

inline long sq(ivec2 a, ivec2 b) {
	long x = (long)a.x - b.x;
	long y = (long)a.y - b.y;
	return x * x + y * y;
}

template<typename RNG>
static ipolygon2 random_polygon(int size, RNG& rng) {
	std::uniform_int_distribution<int> dist(0, 9); //INT_MIN, INT_MAX);
	std::unordered_set<ivec2> points;
	while (points.size() < size) {
		ivec2 p(dist(rng), dist(rng));
		points.insert(p);
	}

	// order points into polygon using traveling salesman heuristic
	ipolygon2 poly(points.begin(), points.end());
	bool done = false;
	while (!done) {
		done = true;
		for (int j : range(size - 2))
			for (int i : range(j + 2, size)) {
				auto a = poly[j];
				auto b = poly[(j + 1) % size];
				auto c = poly[i];
				auto d = poly[(i + 1) % size];
				if (sq(a, b) + sq(c, d) > sq(a, c) + sq(b, d)) {
					std::reverse(poly.begin() + j + 1, poly.begin() + i + 1);
					done = false;
				}
			}
	}
	return poly;
}

inline long edge_area(int a, int b, int c, int d) {
	return (long(a) + b) * (long(c) - d);
}

static long xarea(const ipolygon2& poly) {
	long area = 0;
    auto a = poly.back();
	for (auto b : poly) {
        area += edge_area(a.x, b.x, a.y, b.y);
		a = b;
	}
	return area;
}

static long xarea(ivec2 a, ivec2 b, ivec2 c) {
	long area = 0;
    area += edge_area(a.x, b.x, a.y, b.y);
    area += edge_area(b.x, c.x, b.y, c.y);
    area += edge_area(c.x, a.x, c.y, a.y);
	return area;
}

TEST_CASE("tesselate 1000") {
	std::default_random_engine rng;
	for (auto i : range(3, 1001)) {
		print("Case: {}\n", i);

		ipolygon2 poly;
		long poly_area = 0;
		while (poly_area == 0) {
			poly = random_polygon(i, rng);
			poly_area = xarea(poly);
		}
		auto tess = tesselate(poly);

		print("POLYGON ((");
		for (auto p : poly)
			print("{}, ", p);
		print("{}))\n", poly.front());
	
		print("MULTIPOLYGON ((");
		for (auto t : tess) {
			if (t != tess.front())
				print(", ");
			print("({}, {})", t, t.a);
		}
		print("))\n");
		for (auto t : tess) {
			print("{}\n", (long)xarea(t.a, t.b, t.c));
		}
		
		// verify that sum of areas of all triangles equals polygon area
		long area = 0;
		for (auto m : tess)
			area += std::abs(xarea(m.a, m.b, m.c));
		REQUIRE(std::abs(poly_area) == area);

		// verify that all edges are unique
		std::unordered_set<isegment2> edges;
		for (const itriangle2& triangle : tess)
			for (auto e : triangle.edges()) {
				REQUIRE(edges.count(e) == 0);
				edges.insert(e);
			}

		// all edges should have their opposite except for edges on polygon boundary
		std::unordered_set<isegment2> poly_edges;
		auto a = poly.back();
		for (auto b : poly) {
			poly_edges.insert({a, b});
			a = b;
		}
		for (auto e : edges) {
			auto pc = poly_edges.count(e);
			auto rc = edges.count(e.reversed());
			REQUIRE(pc + rc == 1);
		}
	}
}
