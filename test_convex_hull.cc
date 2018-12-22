#include "convex_hull.h"
#include "is_valid.h"
#include <core/zip.h>
#include "properties.h"
#include "generators.h"
#include <random>
#include <iostream>
#include <algorithm>
#include "mesh_import.h"
#include "is_valid.h"
#include <core/timestamp.h>

using poly_mesh3 = vector<polygon3>;

inline bool operator==(const poly_mesh3& a, const poly_mesh3& b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); i++)
		if (a[i] != b[i])
			return false;
	return true;
}

#include "catch.hpp"

bool less(double4 a, double4 b) {
	if (a.x != b.x)
		return a.x < b.x;
	if (a.y != b.y)
		return a.y < b.y;
	return a.z < b.z;
}

bool less(cspan<double4> v, size_t a, size_t b) {
	for (auto i : range(v.size())) {
		double4 aa = v[(a + i) % v.size()];
		double4 bb = v[(b + i) % v.size()];
		if (!equal(aa, bb))
			return less(aa, bb);
	}
	return false;
};

bool less(cspan<double4> a, cspan<double4> b) {
	if (a.size() != b.size())
		return a.size() < b.size();
	for (auto [aa, bb] : czip(a, b))
		if (!equal(aa, bb))
			return less(aa, bb);
	return false;
};

static aligned_vector<double4> normalize(cspan<double4> v) {
	const auto n = v.size();
	size_t m = 0;
	for (auto i : range<size_t>(1, n))
		if (less(v, i, m))
			m = i;
	aligned_vector<double4> w;
	w.resize(n);
	for (auto i : range(n))
		w[i] = v[(i + m) % n];
	return w;
}

static poly_mesh3 hull(cspan<double4> a) {
	mesh3 m = convex_hull(a);
	// convert mesh3 to ipoly3 with triangles
	poly_mesh3 p;
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

/*static void rotate90_flip_and_shuffle(vector<double4> v) {
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
}*/

TEST_CASE("convex_hull trivial", "[convex_hull]") {
	double4 a = {0, 0, 0, 1};
	double4 b = {1, 0, 0, 1};
	double4 c = {0, 1, 0, 1};
	double4 d = {0, 0, 1, 1};
	REQUIRE(convex_hull({}).empty());
	REQUIRE(convex_hull({a}) == mesh3());
	REQUIRE(convex_hull({a, b}) == mesh3());
	REQUIRE(convex_hull({a, b, c}) == mesh3());

	REQUIRE(convex_hull({a, a, a, a}) == mesh3());
	REQUIRE(convex_hull({a, b, a, b}) == mesh3());
	REQUIRE(convex_hull({a, b, c, a}) == mesh3());
	REQUIRE(convex_hull({a, b, c, b}) == mesh3());
	REQUIRE(convex_hull({a, b, c, c}) == mesh3());
	REQUIRE(convex_hull({a, b, c, d}).size() == 4);
}

TEST_CASE("convex_hull simple", "[convex_hull]") {
	double4 a = {0, 0, 0, 1};
	double4 b = {1000, 0, 0, 1};
	double4 c = {0, 1000, 0, 1};
	double4 d = {0, 0, 1000, 1};

	auto m = hull({a, b, c, d});
	REQUIRE(m.size() == 4);

	REQUIRE(hull({a, b, c, d, double4{1, 1, 1, 1}}) == m);
	REQUIRE(hull({a, b, c, d, double4{0, 1, 1, 1}}) == m);
	REQUIRE(hull({a, b, c, d, double4{1, 0, 1, 1}}) == m);
	REQUIRE(hull({a, b, c, d, double4{1, 1, 0, 1}}) == m);

	REQUIRE(hull({a, b, c, d, double4{-1, 0, 0, 1}}).size() == 4);
}

/*TEST_CASE("convex_hull random points on cube", "[convex_hull]") {
	std::default_random_engine rnd;
	for (int vertices = 4; vertices <= 200; vertices++) {
		aligned_vector<double4> V(vertices);
		for (auto i : range(vertices))
			V[i] = uniform3(rnd, -100, 100);
		mesh3 m = convex_hull(V);
		REQUIRE(IsValid(m) == Validity::OK);
		REQUIRE(is_convex(m));
	}
}*/

TEST_CASE("convex_hull cube", "[convex_hull]") {
	auto m = generate_box(1, 1, 1);
	REQUIRE(is_convex(m));
	REQUIRE(is_aabb(m));
}

class AutoTimestamp {
public:
	AutoTimestamp(string_view name) : _name(name) {}
	~AutoTimestamp() { print("%s took %s ms\n", _begin.elapsed_ms()); }
private:
	string_view _name;
	Timestamp _begin;
};

template<typename Func>
auto measure(string_view name, Func func) {
	AutoTimestamp mm(name);
	return func();
}

#define MEASURE(X) measure(#X, [&](){ return X; } )

template<typename T, typename std::enable_if<std::is_enum<T>::value>::type>
inline void format_e(string& s, string_view spec, T v) {
	format_s(s, "%s(%s)", typeid(T).name(), reinterpret_cast<long>(v));
}

inline void format_e(string& s, string_view spec, Validity v) {
	format_s(s, "Validity(%s)", (int)v);
}

TEST_CASE("convex_hull bunny benchmark", "[!hide][convex_hull]") {
		try{
	std::default_random_engine rnd(0);
	mesh3 mm = MEASURE(load_stl("models/bunny.stl"));
	/*vector<double4> vertices;
	for (const itriangle3& f : mm)
		for (auto i : range(3))
			vertices.push_back(f[i]);
	print("IsValid %s\n", MEASURE(is_valid(mm)));
	mesh3 ch = MEASURE(convex_hull(vertices));
	print("Volume %s\n", MEASURE(volume(ch)));
	print("CenterOfMass %s\n", MEASURE(center_of_mass(ch)));*/
		} catch (char const* e) {
			print("Exception [%s]\n", e);
		}
}
