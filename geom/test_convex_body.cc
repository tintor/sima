#include <catch.hpp>
#include <core/exception.h>
#include <core/timestamp.h>
#include <core/range.h>
#include <core/string_util.h>
#include <geom/convex_body.h>
#include <geom/quaternion.h>

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
	CHECK(expected.size() == actual.size());
	if (expected.size() != actual.size())
		return;

	double dist = std::numeric_limits<double>::infinity();
	for (int shift : range<int>(0, expected.size()))
		for (bool flip : {false,  true})
			minimize(dist, vdist(expected, actual, shift, flip));
	CHECK(dist <= 1e-5);
}

void Check(double3 expected, double3 actual) {
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
	int result = MEASURE(ClassifyConvexConvex(ma, mb, false, normal, contacts, overlap, work1, work2));
	print("result %s, normal %s, contacts %s, overlap %s\n", result, normal, contacts, overlap);

	CHECK(expectedResult == result);
	if (result <= 0 && !equal(expectedNormal, double3{0, 0, 0})) {
		Check(expectedNormal, normal);
	}
	if (result == 0) {
		CheckPoly(expectedContacts, contacts);
	}
	if (result < 0) {
		CHECK(Approx(expectedOverlap).epsilon(1e-6) == overlap);
	}
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
	double3 normal;
	vector<double3> contacts;
	double overlap;

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

// TODO transform contact cases into disjoint by moving away in normal direction, and into penetration by moving closer

// TODO generate two sphere-like meshes and collide them
// TODO generate two capsule-like meshes and collide them

// TODO reuse 2d classify tests for concave polygons

// TODO move shapes into each other until collision

// TODO transformation invariants and swap args invariant

// TODO randomized stress test with two random point clouds
// TODO randomized stress test with face/face contact
// TODO randomized stress test with face/vertex contact
// TODO randomized stress test with edge/edge contact (2 variants)
// TODO randomized stress test with vertex/edge contact
// TODO randomized stress test with vertex/vertex contact
// TODO randomized stress test with face/edge contact

// TODO randomized tests with forced normal along which intervals are in contact, but objects are disjoint!

// TODO benchmarks!
