#include "primitives.hh"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using std::cout;
using std::endl;

TEST_CASE("line_point_squared_distance") {
    REQUIRE(line_point_squared_distance(dvec3(0,0,0), dvec3(9,0,0), dvec3(1,2,0)) == 4);
}

TEST_CASE("segment::nearest") {
    CHECK(segment3(dvec3(0,0,0), dvec3(9,0,0)).nearest(dvec3(1,2,0)) == dvec3(1,0,0));
}

TEST_CASE("Angle") {
	dvec3 a(1, 0, 0);
	dvec3 b(1, 1, 0);
	dvec3 c(cos(M_PI / 6), sin(M_PI / 6), 0);
	dvec3 d(-1, 0, 0);
	CHECK((angle(a, b) - 45_deg) <= 1e-20);
	CHECK((angle(a, c) - 30_deg) <= 1e-20);
	CHECK((angle(a, a) - 0_deg) <= 1e-20);
	CHECK((angle(a, d) - 180_deg) <= 1e-20);
}

TEST_CASE("LineNearest - basic") {
	segment3 p(dvec3(0, 0, 0), dvec3(0, 0, 0));
	segment3 q(dvec3(1, 0, 0), dvec3(1, 0, 0));
	segment3 r(dvec3(1, 0, 0), dvec3(1, 1, 0));

	segment3 pq = segment3::Nearest(p, q).first;
	REQUIRE(pq.a == dvec3(0, 0, 0));
	REQUIRE(pq.b == dvec3(1, 0, 0));

	segment3 pr = segment3::Nearest(p, r).first;
	REQUIRE(pr.a == dvec3(0, 0, 0));
	REQUIRE(pr.b == dvec3(1, 0, 0));

	segment3 rr = segment3::Nearest(r, r).first;
	REQUIRE(rr.a == rr.b);
	REQUIRE(rr.a.x == 1.0);
	REQUIRE(rr.a.y >= 0.0);
	REQUIRE(rr.a.y <= 1.0);
	REQUIRE(rr.a.z == 0.0);
}

segment3 LineNearestNumeric(segment3 p, segment3 q, std::default_random_engine& rnd) {
	std::normal_distribution<real> gauss(0.0, 1.0);
	real snum = 0.5, tnum = 0.5;
	real dnum = l2Norm(p.linear(snum) - q.linear(tnum));
	FOR(j, 1000000) {
		real s = clamp(snum + gauss(rnd));
		real t = clamp(tnum + gauss(rnd));
		const dvec3 A = p.linear(s);
		const dvec3 B = q.linear(t);
		if (l2Norm(A-B) < dnum) {
			dnum = l2Norm(A-B);
			snum = s;
			tnum = t;
		}
	}
	return segment3(p.linear(snum), q.linear(tnum));
}

void TestLineNearest(segment3 p, segment3 q, std::default_random_engine& rnd) {
	segment3 n = segment3::Nearest(p, q).first;
	cout << p.param(n.a) << " " << q.param(n.b) << " computed Angle " << deg(angle(p, n)) << " " << deg(angle(q, n)) << endl;
	real d = l2Norm(n.a - n.b);

	segment3 nn = LineNearestNumeric(p, q, rnd);
	real snum = p.param(nn.a), tnum = q.param(nn.b);
	real dnum = l2Norm(nn.a-nn.b);
	if (dnum < d * 0.99999)
		cout << snum << " " << tnum << " Angle: " << deg(angle(p, nn)) << " " << deg(angle(q, nn)) << endl;
	REQUIRE(dnum >= d * 0.99999);
}

TEST_CASE("LineNearest - random long and short") {
	std::default_random_engine rnd;
	std::uniform_real_distribution<real> uni2(-6, 6);
	FOR(i, 20) {
		cout << "Case " << i << endl;
		dvec3 pa = random_vector(rnd) * 1000.0;
		dvec3 qa = random_vector(rnd) * 1000.0;
		dvec3 pd = random_vector(rnd) * 1000.0 * static_cast<real>(pow(10, uni2(rnd)));
		dvec3 qd = random_vector(rnd) * 1000.0 * static_cast<real>(pow(10, uni2(rnd)));
		segment3 p(pa, pa + pd);
		segment3 q(qa, qa + qd);
		TestLineNearest(p, q, rnd);
	}
}

TEST_CASE("LineNearest - random parallel") {
	std::default_random_engine rnd;
	std::uniform_real_distribution<real> uni2(-1, 1);
	FOR(i, 20) {
		cout << "Case " << i << endl;
		segment3 p(random_vector(rnd) * 1000.0, random_vector(rnd) * 1000.0);
		dvec3 t = random_vector(rnd) * 1000.0;
		dvec3 dir = p.b - p.a;
		segment3 q(p.a + t + uni2(rnd) * dir, p.b + t + uni2(rnd) * dir);
		TestLineNearest(p, q, rnd);
	}
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

TEST_CASE("merge_spheres()") {
	std::default_random_engine rnd;
	real e = 1;
	while (e > 1e-10) {
		sphere a {random_vector(rnd) * 1000.0, 1};
		sphere b {a.center + e * random_direction(rnd), 1};
		sphere c = merge_spheres(a, b);
		REQUIRE(c.radius >= a.radius);
		REQUIRE(c.radius >= b.radius);
		REQUIRE(squared(c.center - a.center) <= squared(c.radius - a.radius) * 1.01);
		REQUIRE(squared(c.center - b.center) <= squared(c.radius - b.radius) * 1.01);
		e *= 0.9;
	}
}
