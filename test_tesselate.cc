#include "tesselate.h"
#include "format.h"
#include "catch.hpp"
#include "range.h"
#include "aabb.h"
#include "string_util.h"
#include <random>
#include "is_valid.h"

template<typename RNG>
static polygon2 random_polygon(int size, RNG& rng) {
	unordered_set<double2, hash_t<double2>, equal_t<double2>> points;
	while (points.size() < size)
		points.insert(uniform2(rng, -1000, 1000));

	// order points into polygon using traveling salesman heuristic
	polygon2 poly(points.begin(), points.end());
	bool done = false;
	while (!done) {
		done = true;
		for (int j : range(1, size))
			for (int i : range(j - 1)) {
				auto c = poly[j];
				auto d = poly[(j + 1) % size];
				auto a = poly[i];
				auto b = poly[(i + 1) % size];
				if (squared(a - b) + squared(c - d) > squared(a - c) + squared(b - d)) {
					// TODO faster if we replace polygon with XOR-linked list
					std::reverse(poly.begin() + i + 1, poly.begin() + j + 1);
					done = false;
				}
			}
	}
	return poly;
}

vector<polygon2> test_cases;
vector<xpolygon2> test_cases2;
mesh2 tess;

#include <fstream>
#include "file.h"
#include "util.h"

struct Setup {
	Setup() {
		std::ifstream is("concave_polygons.txt");
		if (!is.is_open())
			write();
		else {
			is.close();
			read();
		}
	}

	void read() {
		FileReader reader("concave_polygons.txt");
		while (true) {
			string_view line = reader.readline();
			if (line.size() == 0)
				break;
			auto a = split(line);
			polygon2 poly;
			poly.reserve(a.size() / 2);
			for (int i = 0; i < a.size()/2; i++) {
				double x = parse<int>(a[i * 2]);
				double y = parse<int>(a[i * 2 + 1]);
				poly.push_back(double2{x, y});
			}
			xpolygon2 xp;
			xp.add(poly);
			test_cases.push_back(std::move(poly));
			test_cases2.push_back(std::move(xp));
		}

		tess.reserve(test_cases.back().size() + 2);
	}

	void write() {
		std::ofstream os("concave_polygons.txt");
		std::default_random_engine rng;
		for (auto i : range(3, 501)) {
			if (i % 16 == 0)
				print("prepare %s\n", i);
			polygon2 poly;
			while (!IsValid(poly))
				poly = random_polygon(i, rng);
			xpolygon2 xp;
			xp.add(poly);
			test_cases.push_back(poly);
			test_cases2.push_back(std::move(xp));
			for (double2 e : poly)
				os << e.x << ' ' << e.y << ' ';
			os << '\n';
		}

		tess.reserve(test_cases.back().size() + 2);
	}
} setup;

TEST_CASE("tesselate2_500_verify", "[tesselate]") {
	for (const auto& poly : test_cases2) {
		tess.clear();
		tesselate(poly, tess);

		REQUIRE(tess.size() == poly.vertices().size() - 2);

		// verify that sum of areas of all triangles equals polygon area
		double tess_area = 0;
		for (auto m : tess)
			tess_area += area(m);
		double poly_area = area(poly);
		//print("case %d\n", poly.vertices().size());
		REQUIRE(poly_area == tess_area);

		// verify that all edges are unique
		unordered_set<segment2, hash_t<segment2>> edges;
		for (triangle2 triangle : tess)
			for (auto e : Edges(triangle)) {
				REQUIRE(edges.count(e) == 0);
				edges.insert(e);
			}

		// all edges should have their opposite except for edges on polygon boundary
		unordered_set<segment2, hash_t<segment2>> poly_edges;
		for (auto e : Edges(poly))
			poly_edges.insert(e);

		for (auto e : edges) {
			auto pc = poly_edges.count(e);
			auto rc = edges.count(e.reversed());
			REQUIRE(pc + rc == 1);
		}
	}
}

TEST_CASE("tesselate_500_verify", "[tesselate]") {
	for (const auto& poly : test_cases) {
		tess.clear();
		tesselate(poly, tess);

		REQUIRE(tess.size() == poly.size() - 2);

		// verify that sum of areas of all triangles equals polygon area
		double tess_area = 0;
		for (auto m : tess)
			tess_area += area(m);
		double poly_area = area(poly);
		//print("case %d\n", poly.size());
		REQUIRE(poly_area == tess_area);

		// verify that all edges are unique
		unordered_set<segment2, hash_t<segment2>> edges;
		for (triangle2 triangle : tess)
			for (auto e : Edges(triangle)) {
				REQUIRE(edges.count(e) == 0);
				edges.insert(e);
			}

		// all edges should have their opposite except for edges on polygon boundary
		unordered_set<segment2, hash_t<segment2>> poly_edges;
		for (auto e : Edges(poly))
			poly_edges.insert(e);

		for (auto e : edges) {
			auto pc = poly_edges.count(e);
			auto rc = edges.count(e.reversed());
			REQUIRE(pc + rc == 1);
		}
	}
}

TEST_CASE("tesselate2_500_perf", "[!hide][tesselate]") {
	for (const auto& poly : test_cases2) {
		tess.clear();
		tesselate(poly, tess);
	}
}

TEST_CASE("tesselate_500_perf", "[!hide][tesselate]") {
	for (const auto& poly : test_cases) {
		tess.clear();
		tesselate(poly, tess);
	}
}
