#include <catch.hpp>
#include <core/exception.h>
#include <core/timestamp.h>
#include <core/range.h>
#include <core/string_util.h>
#include <geom/convex_body.h>
#include <geom/quaternion.h>

constexpr bool VERBOSE = false;
constexpr int TEST_CASES = 300000;

double3 operator"" _d3(const char* s, size_t n) {
	auto e = split(string_view(s, n), ' ');
	double3 v;
	for (int i : range<int>(3))
		v[i] = parse<double>(e[i]);
	return v;
}

vector<double3> operator"" _vd3(const char* s, size_t n) {
	vector<double3> m;
	for (string_view ss : split(string_view(s, n), ',')) {
		auto e = split(ss, ' ');
		double3 v;
		for (int i : range<int>(3))
			v[i] = parse<double>(e[i]);
		m.push_back(v);
	}
	return m;
}

static double vdist(const vector<double3>& a, const vector<double3>& b, int shift, bool flip) {
	double dist = 0;
	ASSERT_ALWAYS(a.size() == b.size());
	const int N = a.size();
	for (auto i : range(0, N))
		maximize(dist, length(a[((flip ? N - 1 - i : i) + shift) % N] - b[i]));
	return dist;
}

void CheckPoly(const vector<double3>& expected, const vector<double3>& actual) {
	REQUIRE(expected.size() == actual.size());
	if (expected.size() == 0)
		return;

	double dist = std::numeric_limits<double>::infinity();
	for (int shift : range<int>(0, expected.size()))
		for (bool flip : {false,  true})
			minimize(dist, vdist(expected, actual, shift, flip));
	REQUIRE(dist <= 1e-5);
}

void Check(double3 expected, double3 actual) {
	if (VERBOSE)
		print("expected %s, actual %s\n", expected, actual);
	CHECK(length(expected - actual) <= 1e-6);
}

void test(
		const vector<double3>& ca,
		const vector<double3>& cb,
		const int expectedResult,
		const double3 expectedNormal,
		const vector<double3>& expectedContacts,
		const double expectedOverlap) {
	auto ma = GenerateConvexMesh(ca);
	auto mb = GenerateConvexMesh(cb);

	double3 normal = {0, 0, 0};
	double overlap;
	vector<double3> contacts;
	contacts.reserve(8);

	vector<double3> work1;
	work1.reserve(ma.vertices.size() / 2);
	vector<double3> work2;
	work2.reserve(mb.vertices.size() / 2);
	int result = ClassifyConvexConvex(ma, mb, false, normal, contacts, overlap, work1, work2);
	if (VERBOSE)
		print("result %s, normal %s, contacts %s, overlap %s\n", result, normal, contacts, overlap);

	CHECK(expectedResult == result);
	if (result <= 0 && !equal(expectedNormal, double3{0, 0, 0}))
		Check(expectedNormal, normal);
	if (result == 0)
		CheckPoly(expectedContacts, contacts);
	if (result < 0)
		CHECK(Approx(expectedOverlap).epsilon(1e-6) == overlap);
	if (VERBOSE)
		print("\n");
}

vector<double3> translate(const vector<double3>& v, double3 delta) {
	vector<double3> m;
	for (double3 a : v)
		m.push_back(a + delta);
	return m;
}

vector<double3> scale(const vector<double3>& v, double3 scale) {
	vector<double3> m;
	for (double3 a : v)
		m.push_back(a * scale);
	return m;
}

vector<double3> rotate(const vector<double3>& v, double3 axis, double angle) {
	quat q = quat_from_axis_angle(axis, angle);
	vector<double3> m;
	for (double3 a : v)
		m.push_back(quat_rotate(q, a));
	return m;
}

TEST_CASE("convex_body basic", "[convex_body]") {
	vector<double3> ca;
	for (double x : {0, 1})
		for (double y : {0, 1})
			for (double z : {0, 1})
				ca.push_back(double3{x, y, z});

	print("two identical cubes\n");
	test(ca, ca, -1, double3{0, 0, 0}, {}, 1);

	// print("small cube inside big cube\n");
	// test(ca, scale(ca, 0.5), -1, double3{0, 0, 0}, {}, 0.5);

	print("slightly penetrating cubes\n");
	test(ca, translate(ca, double3{0.9, 0.9, 0.95}), -1, double3{0, 0, 1}, {}, 0.05);

	print("vertex / vertex\n");
	test(ca, translate(ca, double3{1, 1, 1}), 0, double3{1, 0, 0}, {double3{1, 1, 1}}, 0);

	print("vertex / edge\n");
	vector<double3> ce = {double3{1, 0.5, 0}};
	for (double y : {0, 1})
		for (double z : {0, 1})
			ce.push_back(double3{2, y, z});
	test(ca, ce, 0, double3{0, 0, 0}, {ce[0]}, 0);

	print("vertex / face\n");
	vector<double3> cd = {double3{1, 0.5, 0.5}};
	for (double y : {0, 1})
		for (double z : {0, 1})
			cd.push_back(double3{2, y, z});
	test(ca, cd, 0, double3{1, 0, 0}, {cd[0]}, 0);

	print("face / face | full face overlap\n");
	test(ca, translate(ca, double3{1, 0, 0}), 0, double3{1, 0, 0},
		"1 0 0, 1 0 1, 1 1 1, 1 1 0"_vd3, 0);

	print("face / face | full face overlap with one cube rotated 45 degrees\n");
	test(ca,
		translate(rotate(translate(ca, double3{1, -0.5, -0.5}), double3{1, 0, 0}, PI/4), double3{0, 0.5, 0.5}),
		0,
		double3{1, 0, 0},
		"1 0 0.292893, 1 0 0.707107, 1 0.292893 1, 1 0.707107 1, "
		"1 1 0.707107, 1 1 0.292893, 1 0.707107 0, 1 0.292893 0"_vd3,
		0);

	print("face / face | small face fully inside larger face\n");
	test(ca, translate(scale(ca, double3{0.5, 0.5, 0.5}), double3{1, 0.25, 0.25}), 0, double3{1, 0, 0},
		"1 0.25 0.25, 1 0.25 0.75, 1 0.75 0.75, 1 0.75 0.25"_vd3, 0);

	print("edge / edge | point in middle\n");
	test(ca, "1 -1 0.5, -1 1 0.5, -1 -1 1, -1 -1 0"_vd3, 0, double3{-0.707107, -0.707107, 0}, "0 0 0.5"_vd3, 0);

	print("edge / edge | full edge overlap\n");
	test(ca, translate(ca, double3{1, 1, 0}), 0, double3{0, 0, 0}, "1 1 0, 1 1 1"_vd3, 0);

	print("edge / edge | smaller edge fully overlaps inside larger edge\n");
	test(ca, translate(scale(ca, "0.5 0.5 0.5"_d3), "1 1 0.25"_d3), 0, double3{0, 0, 0}, "1 1 0.25, 1 1 0.75"_vd3, 0);

	print("edge / edge | edge overlap shifted\n");
	test(ca, translate(ca, double3{1, 1, 0.5}), 0, double3{0, 0, 0}, "1 1 0.5, 1 1 1"_vd3, 0);
}

template<typename T, typename RND>
T uniform_int(RND& rnd, T min, T max) {
	std::uniform_int_distribution dist(min, max);
	return dist(rnd);
}

struct Result {
	int type;
	double3 normal;
	double overlap;
	vector<double3> contacts;
};

template<typename T>
Result Run(const T& ma, const T& mb) {
	vector<double3> work1, work2;
	Result result;
	result.type = ClassifyConvexConvex(ma, mb, false, result.normal, result.contacts, result.overlap, work1, work2);
	return result;
}

Result VerifyInvariants(const vector<double3>& ca, const vector<double3>& cb) {
	auto ma = GenerateConvexMesh(ca);
	auto mb = GenerateConvexMesh(cb);

	auto result = Run(ma, mb);
	if (VERBOSE) {
		print("\nA %s, B %s, result %s, normal %s, contacts %s, overlap %s\n",
			ma.vertices.size(), mb.vertices.size(), result.type, result.normal, result.contacts, result.overlap);
		print("A %s\n", ma.vertices);
		print("B %s\n", mb.vertices);
	}
	REQUIRE(length(result.normal) == Approx(1).margin(1e-4));
	if (result.type == 0)
		REQUIRE(result.contacts.size() > 0);
	// TODO if type == 1 then all vertices of one shape are outside of other shape
	// TODO if any vertex of one shape is inside (or touching) other shape when type == -1
	// TODO check that contacts have no duplicate points!
	// TODO check that all contact points are on boundary of both bodies
	// TODO check that all contact points are in the same plane
	// TODO check that all contact points are forming convex polygon (if >=3 points)
	// TODO check that every contact point is either original OR intersection of two edges

	// TODO if type == 0 then translation of B in normal direction by Tolerance will result in type=1
	// TODO if type == 0 then translation of B in -normal direction by Tolerance will result in type=-1
	// TODO if type == -1 then translation of B in normal direction by overlap will result in type=0

	auto result2 = Run(mb, ma);
	if (VERBOSE) {
		print("\nA %s, B %s, result %s, normal %s, contacts %s, overlap %s\n",
			ma.vertices.size(), mb.vertices.size(), result2.type, result2.normal, result2.contacts, result2.overlap);
		print("C1 %s\n", result.contacts);
		print("C2 %s\n", result2.contacts);
	}
	REQUIRE(result.type == result2.type);
	CheckPoly(result.contacts, result2.contacts);
	return result;
}

template<typename RND>
void Generate(RND& rnd, int pos, int zero, int neg, vector<double3>& ca, int mode = 0) {
	REQUIRE(pos + zero + neg >= 4);
	ca.clear();
	for (auto i : range(0, pos))
		ca.push_back(double3{uniform(rnd, -1, 1), uniform(rnd, -1, 1), uniform(rnd, 1e-3, 1)});
	for (auto i : range(0, zero)) {
		if (mode == 0)
			ca.push_back(double3{uniform(rnd, -1, 1), uniform(rnd, -1, 1), 0});
		if (mode == 1)
			ca.push_back(double3{uniform(rnd, -1, 1), 0, 0});
		if (mode == 2)
			ca.push_back(double3{0, 0, 0});
	}
	for (auto i : range(0, neg))
		ca.push_back(double3{uniform(rnd, -1, 1), uniform(rnd, -1, 1), uniform(rnd, -1, -1e-3)});
	REQUIRE(ca.size() >= 4);
}

TEST_CASE("convex_body random face/face", "[.][convex_body]") {
	std::default_random_engine rnd;
	int count[3] = {0, 0, 0};
	array<int, 20> contacts_count = {0};
	vector<double3> ca, cb;
	for (auto test : range(0, TEST_CASES)) {
		Generate(rnd, uniform_int<int>(rnd, 1, 4), uniform_int<int>(rnd, 3, 5), 0, ca);
		Generate(rnd, 0, uniform_int<int>(rnd, 3, 5), uniform_int<int>(rnd, 1, 4), cb);

		auto result = VerifyInvariants(ca, cb);
		count[result.type + 1] += 1;
		contacts_count[result.contacts.size()] += 1;
		REQUIRE(result.type >= 0);

		// TODO move them closer if disjoint (TODO change "overlap" parameter to distance)
		// TODO move them away if penetrating
	}
	print("disjoint %s, contact %s, overlap %s\n", count[2], count[1], count[0]);
	print("contacts %s\n", contacts_count);
}

TEST_CASE("convex_body random face/edge", "[.][convex_body]") {
	std::default_random_engine rnd;
	int count[3] = {0, 0, 0};
	vector<double3> ca, cb;
	for (auto test : range(0, TEST_CASES)) {
		Generate(rnd, uniform_int<int>(rnd, 1, 4), uniform_int<int>(rnd, 3, 5), 0, ca);
		Generate(rnd, 0, 2, uniform_int<int>(rnd, 2, 5), cb);

		auto result = VerifyInvariants(ca, cb);
		count[result.type + 1] += 1;
		REQUIRE(result.type >= 0);
	}
	print("disjoint %s, contact %s, overlap %s\n", count[2], count[1], count[0]);
}

TEST_CASE("convex_body random face/vertex", "[.][convex_body]") {
	std::default_random_engine rnd;
	int count[3] = {0, 0, 0};
	vector<double3> ca, cb;
	for (auto test : range(0, TEST_CASES)) {
		Generate(rnd, uniform_int<int>(rnd, 1, 4), uniform_int<int>(rnd, 3, 5), 0, ca);
		Generate(rnd, 0, 1, uniform_int<int>(rnd, 3, 5), cb);

		auto result = VerifyInvariants(ca, cb);
		count[result.type + 1] += 1;
		REQUIRE(result.type >= 0);
	}
	print("disjoint %s, contact %s, overlap %s\n", count[2], count[1], count[0]);
}

TEST_CASE("convex_body random edge/edge", "[.][convex_body]") {
	std::default_random_engine rnd;
	int count[3] = {0, 0, 0};
	vector<double3> ca, cb;
	for (auto test : range(0, TEST_CASES)) {
		Generate(rnd, uniform_int<int>(rnd, 2, 5), 2, 0, ca);
		Generate(rnd, 0, 2, uniform_int<int>(rnd, 2, 5), cb);

		auto result = VerifyInvariants(ca, cb);
		count[result.type + 1] += 1;
		REQUIRE(result.type >= 0);
	}
	print("disjoint %s, contact %s, overlap %s\n", count[2], count[1], count[0]);
}

TEST_CASE("convex_body random edge/edge 1d", "[.][convex_body]") {
	std::default_random_engine rnd;
	int count[3] = {0, 0, 0};
	vector<double3> ca, cb;
	for (auto test : range(0, TEST_CASES)) {
		Generate(rnd, uniform_int<int>(rnd, 2, 5), 2, 0, ca, 1);
		Generate(rnd, 0, 2, uniform_int<int>(rnd, 2, 5), cb, 1);

		auto result = VerifyInvariants(ca, cb);
		count[result.type + 1] += 1;
		REQUIRE(result.type >= 0);
	}
	print("disjoint %s, contact %s, overlap %s\n", count[2], count[1], count[0]);
}

TEST_CASE("convex_body random edge/vertex", "[.][convex_body]") {
	std::default_random_engine rnd;
	int count[3] = {0, 0, 0};
	vector<double3> ca, cb;
	for (auto test : range(0, TEST_CASES)) {
		Generate(rnd, uniform_int<int>(rnd, 2, 5), 2, 0, ca, 1);
		Generate(rnd, 0, 1, uniform_int<int>(rnd, 3, 5), cb, 1);

		auto result = VerifyInvariants(ca, cb);
		count[result.type + 1] += 1;
		REQUIRE(result.type >= 0);
	}
	print("disjoint %s, contact %s, overlap %s\n", count[2], count[1], count[0]);
}

TEST_CASE("convex_body random vertex/vertex", "[.][convex_body]") {
	std::default_random_engine rnd;
	int count[3] = {0, 0, 0};
	vector<double3> ca, cb;
	for (auto test : range(0, TEST_CASES)) {
		Generate(rnd, uniform_int<int>(rnd, 3, 3), 1, 0, ca, 2);
		Generate(rnd, 0, 1, uniform_int<int>(rnd, 3, 3), cb, 2);

		auto result = VerifyInvariants(ca, cb);
		count[result.type + 1] += 1;
		REQUIRE(result.type >= 0);
	}
	print("disjoint %s, contact %s, overlap %s\n", count[2], count[1], count[0]);
}

// TODO transform contact cases into disjoint by moving away in normal direction, and into penetration by moving closer

// TODO generate two sphere-like meshes and collide them
// TODO generate two capsule-like meshes and collide them

// TODO reuse 2d classify tests for concave polygons

// TODO move shapes into each other until collision

// TODO transformation invariants and swap args invariant

// TODO randomized stress test with two random point clouds
// TODO randomized stress test with edge/edge contact (2 variants)

// TODO randomized tests with forced normal along which intervals are in contact, but objects are disjoint!

// TODO benchmarks!
