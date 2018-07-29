#include "tesselate.h"
#include "format.h"
#include "catch.hpp"
#include "range.h"
#include "aabb.h"
#include <random>

inline long sq(ivec2 a, ivec2 b) {
	long x = subi(a.x, b.x);
	long y = subi(a.y, b.y);
	return addi(muli(x, x), muli(y, y));
}

template<typename RNG>
static ipolygon2 random_polygon(int size, RNG& rng) {
	std::uniform_int_distribution<int> dist(-1000, 1000); //INT_MIN, INT_MAX);
	std::unordered_set<ivec2, std::hash<ivec2>, equal_t<ivec2>> points;
	while (points.size() < size) {
		ivec2 p{dist(rng), dist(rng)};
		points.insert(p);
	}

	// order points into polygon using traveling salesman heuristic
	ipolygon2 poly(points.begin(), points.end());
	bool done = false;
	while (!done) {
		done = true;
		for (int j : range(1, size))
			for (int i : range(j - 1)) {
				auto c = poly[j];
				auto d = poly[(j + 1) % size];
				auto a = poly[i];
				auto b = poly[(i + 1) % size];
				if (sq(a, b) + sq(c, d) > sq(a, c) + sq(b, d)) {
					// TODO faster if we replace polygon with XOR-linked list
					std::reverse(poly.begin() + i + 1, poly.begin() + j + 1);
					done = false;
				}
			}
	}
	return poly;
}

bool is_valid(const ipolygon2& poly);

std::vector<ipolygon2> test_cases;

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
			std::string_view line = reader.readline();
			if (line.size() == 0)
				break;
			auto a = split(line);
			ipolygon2 poly;
			poly.reserve(a.size() / 2);
			for (int i = 0; i < a.size()/2; i++) {
				int x = parse<int>(a[i * 2]);
				int y = parse<int>(a[i * 2 + 1]);
				poly.push_back(ivec2{x, y});
			}
			test_cases.push_back(std::move(poly));
		}
	}

	void write() {
		std::ofstream os("concave_polygons.txt");
		std::default_random_engine rng;
		for (auto i : range(3, 501)) {
			if (i % 16 == 0)
				print("prepare %s\n", i);
			ipolygon2 poly;
			while (!is_valid(poly))
				poly = random_polygon(i, rng);
			test_cases.push_back(poly);
			for (ivec2 e : poly)
				os << e.x << ' ' << e.y << ' ';
			os << '\n';
		}
	}
} setup;

TEST_CASE("tesselate_500_verify") {
	imesh2 tess;
	tess.reserve(test_cases.back().size() + 2);
	for (const auto& poly : test_cases) {
		tess.clear();
		tesselate(poly, tess);

		REQUIRE(tess.size() == poly.size() - 2);

		// verify that sum of areas of all triangles equals polygon area
		long tess_area = 0;
		for (auto m : tess)
			tess_area += abs(area(m));
		long poly_area = area(poly);
		REQUIRE(std::abs(poly_area) == tess_area);

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
