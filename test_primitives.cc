#include "primitives.h"
#include "catch.hpp"
#include "sphere.h"

TEST_CASE("line_point_squared_distance") {
    REQUIRE(squared_distance(line3(double3{0,0,0}, double3{9,0,0}), double3{1,2,0}) == 4);
}

TEST_CASE("segment::nearest") {
    CHECK(equal(segment3(double3{0,0,0}, double3{9,0,0}).nearest(double3{1,2,0}), double3{1,0,0}));
}

TEST_CASE("Angle") {
	double3 a{1, 0, 0};
	double3 b{1, 1, 0};
	double3 c{cos(M_PI / 6), sin(M_PI / 6), 0};
	double3 d{-1, 0, 0};
	CHECK((angle(a, b) - 45_deg) <= 1e-20);
	CHECK((angle(a, c) - 30_deg) <= 1e-20);
	CHECK((angle(a, a) - 0_deg) <= 1e-20);
	CHECK((angle(a, d) - 180_deg) <= 1e-20);
}

TEST_CASE("LineNearest - basic") {
	segment3 p(double3{0, 0, 0}, double3{0, 0, 0});
	segment3 q(double3{1, 0, 0}, double3{1, 0, 0});
	segment3 r(double3{1, 0, 0}, double3{1, 1, 0});

	segment3 pq = nearest(p, q).first;
	REQUIRE(equal(pq.a, double3{0, 0, 0}));
	REQUIRE(equal(pq.b, double3{1, 0, 0}));

	segment3 pr = nearest(p, r).first;
	REQUIRE(equal(pr.a, double3{0, 0, 0}));
	REQUIRE(equal(pr.b, double3{1, 0, 0}));

	segment3 rr = nearest(r, r).first;
	REQUIRE(equal(rr.a, rr.b));
	REQUIRE(rr.a.x == 1.0);
	REQUIRE(rr.a.y >= 0.0);
	REQUIRE(rr.a.y <= 1.0);
	REQUIRE(rr.a.z == 0.0);
}

segment3 LineNearestNumeric(segment3 p, segment3 q, std::default_random_engine& rnd) {
	std::normal_distribution<double> gauss(0.0, 1.0);
	double snum = 0.5, tnum = 0.5;
	double dnum = length(p.linear(snum) - q.linear(tnum));
	for (auto j : range(1000000)) {
		double s = clamp(snum + gauss(rnd));
		double t = clamp(tnum + gauss(rnd));
		const double3 A = p.linear(s);
		const double3 B = q.linear(t);
		if (length(A - B) < dnum) {
			dnum = length(A-B);
			snum = s;
			tnum = t;
		}
	}
	return segment3(p.linear(snum), q.linear(tnum));
}

void TestLineNearest(segment3 p, segment3 q, std::default_random_engine& rnd) {
	segment3 n = nearest(p, q).first;
	//cout << p.param(n.a) << " " << q.param(n.b) << " computed Angle " << deg(angle(p, n)) << " " << deg(angle(q, n)) << endl;
	double d = length(n.a - n.b);

	segment3 nn = LineNearestNumeric(p, q, rnd);
	double snum = p.param(nn.a), tnum = q.param(nn.b);
	double dnum = length(nn.a-nn.b);
	//if (dnum < d * 0.99999)
	//	cout << snum << " " << tnum << " Angle: " << deg(angle(p, nn)) << " " << deg(angle(q, nn)) << endl;
	REQUIRE(dnum >= d * 0.99999);
}

TEST_CASE("LineNearest - random long and short") {
	std::default_random_engine rnd;
	std::uniform_real_distribution<double> uni2(-6, 6);
	for (auto i : range(20)) {
		//cout << "Case " << i << endl;
		double3 pa = uniform3(rnd, -1, 1) * 1000.0;
		double3 qa = uniform3(rnd, -1, 1) * 1000.0;
		double3 pd = uniform3(rnd, -1, 1) * 1000.0 * static_cast<double>(pow(10, uni2(rnd)));
		double3 qd = uniform3(rnd, -1, 1) * 1000.0 * static_cast<double>(pow(10, uni2(rnd)));
		segment3 p(pa, pa + pd);
		segment3 q(qa, qa + qd);
		TestLineNearest(p, q, rnd);
	}
}

TEST_CASE("LineNearest - random parallel") {
	std::default_random_engine rnd;
	std::uniform_real_distribution<double> uni2(-1, 1);
	for (auto i : range(20)) {
		//cout << "Case " << i << endl;
		segment3 p(uniform3(rnd, -1, 1) * 1000.0, uniform3(rnd, -1, 1) * 1000.0);
		double3 t = uniform3(rnd, -1, 1) * 1000.0;
		double3 dir = p.b - p.a;
		segment3 q(p.a + t + uni2(rnd) * dir, p.b + t + uni2(rnd) * dir);
		TestLineNearest(p, q, rnd);
	}
}

double3 vec3(double x, double y, double z) {
	return {x, y, z};
}

TEST_CASE("distance(triangle3, triangle3)") {
    // Planar cases:
    // vertex-vertex
    CHECK(distance(triangle3(vec3(0,0,0), vec3(2,0,0), vec3(0,2,0)), triangle3(vec3(5,0,0), vec3(7,0,0), vec3(5,2,0))) == 3);
    // common edge
    CHECK(distance(triangle3(vec3(0,0,0), vec3(2,0,0), vec3(0,2,0)), triangle3(vec3(0,0,0), vec3(2,0,0), vec3(0,-2,0))) == 0);

    // Parallel cases:
    // face-face overlap
    CHECK(distance(triangle3(vec3(0,0,0), vec3(2,0,0), vec3(0,2,0)), triangle3(vec3(0,0,1), vec3(2,0,1), vec3(0,2,1))) == 1);
}

