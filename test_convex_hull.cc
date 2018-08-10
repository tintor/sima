#include "convex_hull.h"
#include "is_valid.h"
#include "zip.h"
#include "properties.h"
#include "generators.h"
#include <random>
#include <iostream>
#include <algorithm>
#include "mesh_import.h"
#include "is_valid.h"
#include "timestamp.h"

using ipoly_mesh3 = std::vector<ipolygon3>;

inline bool operator==(const ipoly_mesh3& a, const ipoly_mesh3& b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); i++)
		if (a[i] != b[i])
			return false;
	return true;
}

#include "catch.hpp"

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
		if (!equal(aa, bb))
			return less(aa, bb);
	}
	return false;
};

bool less(const std::vector<ivec3>& a, const std::vector<ivec3>& b) {
	if (a.size() != b.size())
		return a.size() < b.size();
	for (auto [aa, bb] : czip(a, b))
		if (!equal(aa, bb))
			return less(aa, bb);
	return false;
};

static std::vector<ivec3> normalize(const std::vector<ivec3>& v) {
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

static ipoly_mesh3 hull(std::vector<ivec3> a) {
	imesh3 m = convex_hull(a);
	// convert imesh3 to ipoly3 with triangles
	ipoly_mesh3 p;
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

static void rotate90_flip_and_shuffle(std::vector<ivec3> v) {
	std::default_random_engine rnd;
	auto m = hull(v);
	for (auto i : range(100)) {
		std::shuffle(v.begin(), v.end(), rnd);
		// TODO apply random (negate any axis / swap any two coordinates)
		auto p = hull(v);
		// TODO make sure to do reverse transform on p
		REQUIRE(m.size() == p.size());
		for (size_t i = 0; i < m.size(); i++)
			REQUIRE(m[i] == p[i]);
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

	auto m = hull({a, b, c, d});
	REQUIRE(m.size() == 4);

	REQUIRE(hull({a, b, c, d, ivec3{1, 1, 1}}) == m);
	REQUIRE(hull({a, b, c, d, ivec3{0, 1, 1}}) == m);
	REQUIRE(hull({a, b, c, d, ivec3{1, 0, 1}}) == m);
	REQUIRE(hull({a, b, c, d, ivec3{1, 1, 0}}) == m);

	REQUIRE(hull({a, b, c, d, ivec3{-1, 0, 0}}).size() == 4);
}

// uniform inside a cube
template<typename RandomEngine>
ivec3 random_vector(RandomEngine& rnd, int a, int b) {
	std::uniform_int_distribution<int> dist(a, b);
	return ivec3{dist(rnd), dist(rnd), dist(rnd)};
}

TEST_CASE("convex_hull random points on cube") {
	std::default_random_engine rnd;
	for (int vertices = 4; vertices <= 200; vertices++) {
		std::vector<ivec3> V(vertices);
		for (auto i : range(vertices))
			V[i] = random_vector(rnd, -100, 100);
		imesh3 m = convex_hull(V);
		REQUIRE(is_valid(m) == Validity::OK);
		REQUIRE(is_convex(m));
	}
}

TEST_CASE("convex_hull cube") {
	auto m = generate_box(1, 1, 1);
	REQUIRE(is_convex(m));
	REQUIRE(is_aabb(m));
}

class AutoTimestamp {
public:
	AutoTimestamp(std::string_view name) : _name(name) {}
	~AutoTimestamp() { print("%s took %s ms\n", _begin.elapsed_ms()); }
private:
	std::string_view _name;
	Timestamp _begin;
};

template<typename Func>
auto measure(std::string_view name, Func func) {
	AutoTimestamp mm(name);
	return func();
}

#define MEASURE(X) measure(#X, [&](){ return X; } )

template<typename T, typename std::enable_if<std::is_enum<T>::value>::type>
inline void format_e(std::string& s, std::string_view spec, T v) {
	format_s(s, "%s(%s)", typeid(T).name(), reinterpret_cast<long>(v));
}

inline void format_e(std::string& s, std::string_view spec, Validity v) {
	format_s(s, "Validity(%s)", (int)v);
}

TEST_CASE("convex_hull bunny benchmark", "[!hide]") {
		try{
	std::default_random_engine rnd(0);
	imesh3 mm = MEASURE(load_stl("models/bunny.stl", 1000));
	/*std::vector<ivec3> vertices;
	for (const itriangle3& f : mm)
		for (auto i : range(3))
			vertices.push_back(f[i]);
	print("IsValid %s\n", MEASURE(is_valid(mm)));
	imesh3 ch = MEASURE(convex_hull(vertices));
	print("Volume %s\n", MEASURE(volume(ch)));
	print("CenterOfMass %s\n", MEASURE(center_of_mass(ch)));*/
		} catch (char const* e) {
			print("Exception [%s]\n", e);
		}
}
