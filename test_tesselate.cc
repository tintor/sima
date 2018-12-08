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
mesh2 tess;

#include <fstream>
#include "file.h"
#include "util.h"

string poly_to_str(span<const double2> p) {
	string s;
	for (double2 e : p) {
		format(s, "%f", e.x);
		s += ' ';
		format(s, "%f", e.y);
		s += ' ';
	}
	return s;
}

polygon2 str_to_poly(string_view s) {
	auto a = split(s);
	polygon2 p;
	p.reserve(a.size() / 2);
	for (int i = 0; i < a.size()/2; i++) {
		double x = parse<double>(a[i * 2]);
		double y = parse<double>(a[i * 2 + 1]);
		p.push_back(double2{x, y});
	}
	return p;
}

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
			polygon2 poly = str_to_poly(line);
			if (!IsValid(poly))
				THROW(runtime_error, "polygon not valid");
			test_cases.push_back(std::move(poly));
		}

		tess.reserve(test_cases.back().size() + 2);
	}

	void write() {
		std::ofstream os("concave_polygons.txt");
		std::default_random_engine rng;
		for (auto i : range(3, 701)) {
			if (i % 16 == 0)
				print("prepare %s\n", i);
			polygon2 poly;
			while (!IsValid(poly))
				poly = random_polygon(i, rng);
			test_cases.push_back(poly);
			os << poly_to_str(poly) << '\n';
		}

		tess.reserve(test_cases.back().size() + 2);
	}
} setup;

TEST_CASE("tesselateSimple_verify", "[tesselate]") {
	for (const auto& poly : test_cases) {
		tess.clear();
		tesselateSimple(poly, tess);

		REQUIRE(tess.size() == poly.size() - 2);

		// verify that sum of areas of all triangles equals polygon area
		double tess_area = 0;
		for (auto m : tess)
			tess_area += area(m);
		double poly_area = area(poly);
		print("case %d\n", poly.size());
		REQUIRE(poly_area == Approx(tess_area).margin(1e-6));

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

/*TEST_CASE("tesselate_500_verify", "[tesselate]") {
	for (const auto& poly : test_cases) {
		tess.clear();
		tesselate(poly, tess);

		REQUIRE(tess.size() == poly.size() - 2);

		// verify that sum of areas of all triangles equals polygon area
		double tess_area = 0;
		for (auto m : tess)
			tess_area += area(m);
		double poly_area = area(poly);
		print("case %d\n", poly.size());
		REQUIRE(poly_area == Approx(tess_area).margin(1e-6));

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
}*/

extern long t_deck;
extern long t_pip;
extern long t_relate;
extern long t_erase;

TEST_CASE("tesselateSimple_perf", "[!hide][tesselate]") {
	t_deck = t_pip = t_relate = t_erase = 0;
	for (const auto& poly : test_cases) {
		tess.clear();
		tesselateSimple(poly, tess);
	}
	print("t_deck %d\n", t_deck / 1000000);
	print("t_pip %d\n", t_pip / 1000000);
	print("t_relate %d\n", t_relate / 1000000);
	print("t_erase %d\n", t_erase / 1000000);
}

/*TEST_CASE("tesselate_500_perf", "[!hide][tesselate]") {
	for (const auto& poly : test_cases) {
		tess.clear();
		tesselate(poly, tess);
	}
}*/
