#include "convex_hull.h"
#include "zip.h"
#include "catch.hpp"
#include <random>
#include <iostream>
#include <algorithm>

using ipoly3 = std::vector<std::vector<ivec3>>;

bool less(ivec3 a, ivec3 b) {
	if (a.x != b.x)
		return a.x < b.x;
	if (a.y != b.y)
		return a.y < b.y;
	return a.z < b.z;
}

bool less(const std::vector<ivec3>& v, size_t a, size_t b) {
	for (auto i : range(v.size())) {
		ivec3 aa = v[(a + i) % v.size()];
		ivec3 bb = v[(b + i) % v.size()];
		if (aa != bb)
			return less(aa, bb);
	}
	return false;
};

bool less(const std::vector<ivec3>& a, const std::vector<ivec3>& b) {
	if (a.size() != b.size())
		return a.size() < b.size();
	for (auto [aa, bb] : czip(a, b))
		if (aa != bb)
			return less(aa, bb);
	return false;
};

std::vector<ivec3> normalize(const std::vector<ivec3>& v) {
	const auto n = v.size();
	size_t m = 0;
	for (auto i : range<size_t>(1, n))
		if (less(v, i, m))
			m = i;
	std::vector<ivec3> w;
	w.resize(n);
	for (auto i : range(n))
		w[i] = v[(i + m) % n];
	return w;
}

ipoly3 hull(std::vector<ivec3> a) {
	imesh3 m = convex_hull(a);
	// convert imesh3 to ipoly3 with triangles
	ipoly3 p;
	for (auto f : m)
		p.push_back({ f.a, f.b, f.c });
	// merge touching coplanar faces in ipoly3
	
		
	// sort vertices in each face in ipoly3
	for (auto& f : p)
		normalize(f);
	// sort faces in ipoly3
	// TODO std::sort(p.begin(), p.end(), less);
	return p;	
}

void rotate90_flip_and_shuffle(std::vector<ivec3> v) {
	std::default_random_engine rnd;
	ipoly3 m = hull(v);
	for (auto i : range(100)) {	
		std::shuffle(v.begin(), v.end(), rnd);
		// TODO apply random (negate any axis / swap any two coordinates)
		ipoly3 p = hull(v);
		// TODO make sure to do reverse transform on p
		REQUIRE(m == p);
	}
}

TEST_CASE("convex_hull trivial") {
	ivec3 a = {0, 0, 0};
	ivec3 b = {1, 0, 0};
	ivec3 c = {0, 1, 0};
	ivec3 d = {0, 0, 1};
	REQUIRE(convex_hull(std::vector<ivec3>{}) == imesh3());
	REQUIRE(convex_hull(std::vector{a}) == imesh3());
	REQUIRE(convex_hull(std::vector{a, b}) == imesh3());
	REQUIRE(convex_hull(std::vector{a, b, c}) == imesh3());

	REQUIRE(convex_hull(std::vector{a, a, a, a}) == imesh3());
	REQUIRE(convex_hull(std::vector{a, b, a, b}) == imesh3());
	REQUIRE(convex_hull(std::vector{a, b, c, a}) == imesh3());
	REQUIRE(convex_hull(std::vector{a, b, c, b}) == imesh3());
	REQUIRE(convex_hull(std::vector{a, b, c, c}) == imesh3());
	REQUIRE(convex_hull(std::vector{a, b, c, d}).size() == 4);
}


TEST_CASE("convex_hull simple") {
	ivec3 a = {0, 0, 0};
	ivec3 b = {1000, 0, 0};
	ivec3 c = {0, 1000, 0};
	ivec3 d = {0, 0, 1000};

	ipoly3 m = hull({a, b, c, d});
	REQUIRE(m.size() == 4);

	REQUIRE(hull({a, b, c, d, ivec3(1, 1, 1)}) == m);
	REQUIRE(hull({a, b, c, d, ivec3(0, 1, 1)}) == m);
	REQUIRE(hull({a, b, c, d, ivec3(1, 0, 1)}) == m);
	REQUIRE(hull({a, b, c, d, ivec3(1, 1, 0)}) == m);

	// TODO REQUIRE(hull({a, b, c, d, ivec3(-1, 0, 0)}).size() == 6);
}

// uniform inside a cube
template<typename RandomEngine>
ivec3 random_vector(RandomEngine& rnd, int a = -100, int b = 100) {
	std::uniform_int_distribution<int> dist(a, b);
	return ivec3(dist(rnd), dist(rnd), dist(rnd));
}

/*TEST_CASE("convex_hull random points on cube") {
	std::default_random_engine rnd;
	for (int vertices = 4; vertices <= 40; vertices++) {
		std::vector<ivec3> V(vertices);
		for (auto i : range(vertices))
			V[i] = random_vector(rnd);
		imesh3 m = convex_hull(V);
		REQUIRE(is_convex(m));
		std::cout << vertices << " " << m.size() << std::endl;
	}
}*/

// uniform on the unit sphere surface
template<typename RandomEngine>
ivec3 random_direction(RandomEngine& rnd) {
	std::normal_distribution<double> dist(0.0, 1.0);
	dvec3 v = 1000000.0 * normalize(dvec3(dist(rnd), dist(rnd), dist(rnd)));
	return { std::round(v.x), std::round(v.y), std::round(v.z) };
}

double volume(const imesh3& mesh) {
	double volume = 0;
	for (auto m : mesh) {
		double s = 0, z = 0;
		for (auto [a, b] : m.edges()) {
			z += b.z;
			s += ((double)a.y + b.y) * ((double)a.x - b.x);
		}
		volume += z * s;
	}
	return volume / 6;
}

/*TEST_CASE("convex_hull random points on sphere") {
	std::default_random_engine rnd;
	int vertices = 1000;
	std::vector<ivec3> V(vertices);
	for (auto i : range(vertices))
		V[i] = random_direction(rnd);
	imesh3 m = convex_hull(V);
	double v = 4.0 / 3 * M_PI * 1e18;
	REQUIRE(m.size() == vertices * 2 - 4);
	REQUIRE(abs(volume(m) - v) / v < 0.15);
}*/
