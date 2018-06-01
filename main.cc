#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <functional>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "tinyformat.h"
#include "auto.h"

using std::pair;
using std::cerr;
using std::cout;
using std::endl;
using std::size_t;
using std::min;
using std::max;

using namespace tfm;

template<typename T> constexpr
T min(T a, T b, T c) { return min(a, min(b, c)); }
template<typename T> constexpr
T max(T a, T b, T c) { return max(a, max(b, c)); }
template<typename T> constexpr
T min(T a, T b, T c, T d) { return min(min(a, b), min(c, d)); }
template<typename T> constexpr
T max(T a, T b, T c, T d) { return max(max(a, b), max(c, d)); }

#define FOR(I, N) for (decltype(N) I = 0; I < N; I++)
#define FOR_EACH(A, B) for (auto& A : B)
#define FOR_EACH_EDGE(A, B, V) for (auto *B = V.begin(), *A = V.begin() + 2; B < V.begin() + 3; A = B++)

using real = double;

#define GLM_ENABLE_EXPERIMENTAL
#define GLM_FORCE_CXX14
#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtx/quaternion.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_access.hpp"
#include "glm/gtx/norm.hpp"
using namespace glm;

static_assert(sizeof(dvec3) == sizeof(real) * 3);

constexpr real squared(real a) { return a * a; }
real squared(dvec3 a) { return dot(a, a); }
real squared(vec4 a) { return dot(a, a); }

constexpr long double operator "" _deg(long double a) { return a * (M_PI / 180); }
constexpr long double operator "" _deg(unsigned long long a) { return a * (M_PI / 180); }

std::string deg(long double r) { return format("%llg deg", r * (180 / M_PI)); }

auto angle(dvec3 a, dvec3 b) {
	// atan2 is numericaly better than acos when angle is very small

	// make the vectors equal length first (without division)
	real aa = l2Norm(a), bb = l2Norm(b);
	a *= bb;
	b *= aa;
	return 2 * atan2(l2Norm(a - b), l2Norm(a + b));
}

constexpr real clamp(real t, real min = 0, real max = 1) {
	if (t < min) return min;
	if (t > max) return max;
	return t;
}

struct ray3
{
    // TODO is dir constrained to UnitVector?
    const dvec3 origin, dir;
};

// infinite line vs. point
real line_point_squared_distance(dvec3 a, dvec3 b, dvec3 p) {
	dvec3 ba = b - a, pa = p - a;
    real t = dot(pa, ba) / dot(ba, ba);
    return squared(pa - ba * t);
}

TEST_CASE("line_point_squared_distance") {
    REQUIRE(line_point_squared_distance(dvec3(0,0,0), dvec3(9,0,0), dvec3(1,2,0)) == 4);
}

// line segment
struct segment3 {
	dvec3 a, b;

	segment3() { }
	constexpr segment3(dvec3 a, dvec3 b) : a(a), b(b) { }

	bool operator==(segment3 p) const { return a == p.a && b == p.b; }
	bool operator!=(segment3 p) const { return a != p.a || b != p.b; }

	constexpr segment3 reverse() const { return segment3(b, a); }
	dvec3 linear(real t) const { return a + (b - a) * t; }

	// Inverse of linear
	real param(dvec3 p) const {
		dvec3 d = b - a;
		return dot(p - a, d) / dot(d, d);
	}

	// Return point on line segment that is nearest to P
	dvec3 nearest(dvec3 p) const {
		dvec3 d = b - a;
		real t = dot(p - a, d);
		if (t <= 0)
			return a;
		if (t >= dot(d, d))
			return b;
		return a + d * (t / dot(d, d));
	}

	// Tested and stable!
	enum class NearestCase {
		RayRay = 0,
		RayPoint = 1,
		PointPoint = 2,
	};
	static pair<segment3, NearestCase> Nearest(segment3 p, segment3 q) {
		dvec3 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
		real aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
		constexpr real inf = std::numeric_limits<real>::max();
		constexpr real tiny = 1e-8;

		// ray/ray
		real d = aa * bb - ab * ab;
		real s = ab * bc - bb * ac;
		real t = aa * bc - ab * ac;
		// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
		if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
				|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0))
			return pair<segment3, NearestCase>(segment3(p.a + A * (s / d), q.a + B * (t / d)), NearestCase::RayRay);

		// ray/endpoint
		real s0 = (aa >= tiny) ? -ac / aa : -1;
		real s1 = (aa >= tiny) ? (ab - ac) / aa : -1;
		real t0 = (bb >= tiny) ? bc / bb : -1;
		real t1 = (bb >= tiny) ? (ab + bc) / bb : -1;

		real d1 = (0 <= s0 && s0 <= 1) ? squared(C + A*s0) : inf;
		real d2 = (0 <= s1 && s1 <= 1) ? squared(C - B + A*s1) : inf;
		real d3 = (0 <= t0 && t0 <= 1) ? squared(B*t0 - C) : inf;
		real d4 = (0 <= t1 && t1 <= 1) ? squared(B*t1 - C - A) : inf;

		// endpoint/endpoint
		real d5 = squared(C);
		real d6 = squared(C + A);
		real d7 = squared(C - B);
		real d8 = squared(C + A - B);

		real dm = std::min(min(d1, d2, d3, d4), min(d5, d6, d7, d8));

		if (d1 == dm)
			return pair(segment3(p.a + A*s0, q.a), NearestCase::RayPoint);
		if (d2 == dm)
			return pair(segment3(p.a + A*s1, q.b), NearestCase::RayPoint);
		if (d3 == dm)
			return pair(segment3(p.a, q.a + B*t0), NearestCase::RayPoint);
		if (d4 == dm)
			return pair(segment3(p.b, q.a + B*t1), NearestCase::RayPoint);
		if (d5 == dm)
			return pair(segment3(p.a, q.a), NearestCase::PointPoint);
		if (d6 == dm)
			return pair(segment3(p.b, q.a), NearestCase::PointPoint);
		if (d7 == dm)
			return pair(segment3(p.a, q.b), NearestCase::PointPoint);
		if (d8 == dm)
			return pair(segment3(p.b, q.b), NearestCase::PointPoint);

		throw new std::exception();
	}

	static real squared_distance(segment3 p, segment3 q) {
		dvec3 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
		real aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
		constexpr real inf = std::numeric_limits<real>::max();
		constexpr real tiny = 1e-8;

		// ray/ray
		real d = aa * bb - ab * ab;
		real s = ab * bc - bb * ac;
		real t = aa * bc - ab * ac;
		// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
		if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
				|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0))
			return squared(C + A * (s / d) - B * (t / d));

		// ray/endpoint
		real s0 = (aa >= tiny) ? -ac / aa : -1;
		real s1 = (aa >= tiny) ? (ab - ac) / aa : -1;
		real t0 = (bb >= tiny) ? bc / bb : -1;
		real t1 = (bb >= tiny) ? (ab + bc) / bb : -1;

		real d1 = (0 <= s0 && s0 <= 1) ? squared(C + A*s0) : inf;
		real d2 = (0 <= s1 && s1 <= 1) ? squared(C - B + A*s1) : inf;
		real d3 = (0 <= t0 && t0 <= 1) ? squared(B*t0 - C) : inf;
		real d4 = (0 <= t1 && t1 <= 1) ? squared(B*t1 - C - A) : inf;

		// endpoint/endpoint
		real d5 = squared(C);
		real d6 = squared(C + A);
		real d7 = squared(C - B);
		real d8 = squared(C + A - B);

		return std::min(min(d1, d2, d3, d4), min(d5, d6, d7, d8));
	}
};

TEST_CASE("segment::nearest") {
    CHECK(segment3(dvec3(0,0,0), dvec3(9,0,0)).nearest(dvec3(1,2,0)) == dvec3(1,0,0));
}

auto angle(segment3 p, segment3 q) {
	return angle(p.b - p.a, q.b - q.a);
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

dvec3 UniVec(std::default_random_engine& rnd) {
	std::uniform_real_distribution<real> uni(-1000, 1000);
	return dvec3(uni(rnd), uni(rnd), uni(rnd));
}

TEST_CASE("LineNearest - random long and short") {
	std::default_random_engine rnd;
	std::uniform_real_distribution<real> uni2(-6, 6);
	FOR(i, 20) {
		cout << "Case " << i << endl;
		dvec3 pa = UniVec(rnd);
		dvec3 qa = UniVec(rnd);
		dvec3 pd = UniVec(rnd) * static_cast<real>(pow(10, uni2(rnd)));
		dvec3 qd = UniVec(rnd) * static_cast<real>(pow(10, uni2(rnd)));
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
		segment3 p(UniVec(rnd), UniVec(rnd));
		dvec3 t = UniVec(rnd);
		dvec3 dir = p.b - p.a;
		segment3 q(p.a + t + uni2(rnd) * dir, p.b + t + uni2(rnd) * dir);
		TestLineNearest(p, q, rnd);
	}
}

struct triangle3 {
	dvec3 a, b, c;

    triangle3() { }
	constexpr triangle3(dvec3 a, dvec3 b, dvec3 c) : a(a), b(b), c(c) { }

    constexpr const dvec3* begin() const { return &a; }
    constexpr dvec3 operator[](int i) const { return (&a)[i]; }

	constexpr void operator-=(dvec3 v) {
		a -= v;
		b -= v;
		c -= v;
	}
};

typedef std::vector<triangle3> Mesh3d;

dvec3 Normal(const dvec3& a, const dvec3& b, const dvec3& c) {
	// TODO make sure it works when one edge is very small compared to other two
	// TODO make sure it works when one edge is almost equal to the sum of other two
	return glm::cross(b - a, c - a);
}

dvec3 Normal(const triangle3& p) {
	return Normal(p.a, p.b, p.c);
}

// How close a point needs to be to a plane to be considered on the plane?
constexpr real PlanarEpsilon = 1e-6;

// How close two features need to be to be considered touching?
// TODO can it be same as PlanarEpsilon?
constexpr real ContactEpsilon = 10 * PlanarEpsilon;

struct plane {
	const dvec3 normal; // unit vector
	const real d;

	plane(const dvec3& a, const dvec3& b, const dvec3& c)
		: normal(normalize(Normal(a, b, c))), d(dot(-normal, a)) { }

	plane(const triangle3& m) : plane(m.a, m.b, m.c) {}

	// Signed distance!
	real distance(const dvec3& p) const {
		return dot(normal, p) + d;
	}

	// > 0, if P is on the positive side of plane ABC (right hand rule)
	static real sign(const dvec3& a, const dvec3& b, const dvec3& c, const dvec3& p) {
		return dot(Normal(a, b, c), p - a);
	}
};

constexpr bool intersects(const segment3& q, const plane& p) {
	real a = p.distance(q.a);
	real b = p.distance(q.b);
	real e = PlanarEpsilon;
	return (a <= e || b <= e) && (a >= -e || b >= -e);
}

constexpr bool intersects(const triangle3& q, const plane& p) {
	real a = p.distance(q.a);
	real b = p.distance(q.b);
	real c = p.distance(q.c);
	real e = PlanarEpsilon;
	return (a <= e || b <= e || c <= e) && (a >= -e || b >= -e || c >= -e);
}

struct plucker {
	dvec3 u, v;

	static plucker Line(const dvec3& a, const dvec3& b) {
		return plucker{b - a, cross(b, a)};
	}

	// <0 Clockwise (if you look in direction of one line, other will go CW around it)
 	// =0 Intersect or Parallel
  	// >0 Counterclockwise
	real crossing(plucker p) const {
		return dot(u, p.v) + dot(v, p.u);
	}
};

dvec3 mini(const triangle3& p) {
	dvec3 e;
	FOR(i, 3)
		e[i] = min(p.a[i], p.b[i], p.c[i]);
	return e;
}

dvec3 maxi(const triangle3& p) {
	dvec3 e;
	FOR(i, 3)
		e[i] = max(p.a[i], p.b[i], p.c[i]);
	return e;
}

bool DisjointIntervals(real amin, real amax, real bmin, real bmax) {
	return amax + PlanarEpsilon < bmin || bmax + PlanarEpsilon < amin;
}

// Valid if triangles are not intersecting, except in one shared edge or one shared vertex
bool AreValidMeshFaces(const triangle3& a, const triangle3& b) {
	// Axis check for early exit
	dvec3 amin = mini(a), amax = maxi(a);
	dvec3 bmin = mini(b), bmax = maxi(b);
	FOR(i, 3)
		if (DisjointIntervals(amin[i], amax[i], bmin[i], bmax[i]))
			return true;

	// Plane check
	if (!intersects(a, plane(b)) || !intersects(b, plane(a)))
		return true;

	// Non-planar case
	int match[3] = {-1, -1, -1};
	FOR(i, 3) FOR(j, 3)
		if (a[i] == b[j])
			match[i] = j;
	int count = 0;
	FOR(i, 3)
		if (match[i] >= 0)
			count += 1;
	if (count == 1) {

	}

	// Planar case
	// TODO there is some intersection:
	// OK case is if intersection is one vertex of both A and B (and no overlap)
	// OK case is if intersection is one edge of both A and B (and no overlap)
	return true;
}

class UnionFind {
public:
	UnionFind() : m_parent(this), m_rank(0) { }

	void Union(UnionFind& b) {
		UnionFind* pa = Find();
		UnionFind* pb = b.Find();

		if (pa->m_rank < pb->m_rank) {
			pa->m_parent = pb;
		} else if (pa->m_rank > pb->m_rank) {
			pb->m_parent = pa;
		} else {
			pa->m_rank += 1;
			pb->m_parent = pa;
		}
	}

	UnionFind* Find() {
		// Path halving faster than path compression (from Wikipedia)
		UnionFind* x = this;
		while (x->m_parent != x) {
			x->m_parent = x->m_parent->m_parent;
			x = x->m_parent;
		}
  		return x;
	}

private:
	UnionFind* m_parent;
	int m_rank;
};

namespace std {
	// Hash function for Eigen matrix and vector.
	// The code is from `hash_combine` function of the Boost library. See
	// http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine .
	template<> struct hash<dvec3> {
		size_t operator()(const dvec3& a) const {
			size_t seed = 0;
			FOR(i, 3) seed ^= std::hash<real>()(a[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			return seed;
		}
    };

  template <typename A, typename B> struct hash<pair<A, B>> {
    size_t operator()(const pair<A, B>& x) const {
		size_t seed = 0;
		seed ^= std::hash<A>()(x.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= std::hash<B>()(x.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	    return seed;
    }
  };

  template<> struct hash<segment3> {
    size_t operator()(const segment3& x) const {
		size_t seed = 0;
		seed ^= std::hash<dvec3>()(x.a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= std::hash<dvec3>()(x.b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	    return seed;
    }
  };
}

enum class Validity {
	OK = 0,
	TooFewFaces = 1,
	EdgeTooShort = 2,
	OpenEdge = 3,
	SeparateComponents = 4,
	SelfIntersection = 5,
};

Validity IsValid(const Mesh3d& mesh) {
	if (mesh.size() < 4)
		return Validity::TooFewFaces;

	// Minimal length of any edge is 10xPlanarEpsilon
	FOR(i, mesh.size())
		FOR_EACH_EDGE(a, b, mesh[i])
			if (squared(*a - *b) <= 100 * squared(PlanarEpsilon))
				return Validity::EdgeTooShort;

	// TODO any two different vertices can't be closer than PlanarEpsilon
	// TODO any vertex can't be closer than PlanarEpsilon to any other edge

	// Extract all edges
	std::unordered_map<segment3, int> edge_to_face;
	FOR(i, mesh.size())
		FOR_EACH_EDGE(p, q, mesh[i])
 			edge_to_face[segment3(*p, *q)] = i;

	// Every edge must appear exactly twice (in opposite directions)
	FOR_EACH(e, edge_to_face)
		if (edge_to_face.count(e.first.reverse()) == 0)
			return Validity::OpenEdge;

	// All triangles must form a single component
	std::vector<UnionFind> component(mesh.size());
	FOR(i, mesh.size())
		FOR_EACH_EDGE(p, q, mesh[i])
			component[edge_to_face[segment3(*q, *p)]].Union(component[i]);
	FOR(i, component.size())
		if (component[i].Find() != component[0].Find())
			return Validity::SeparateComponents;

	// No self intersections
	FOR(i, mesh.size())
		FOR(j, i)
			if (!AreValidMeshFaces(mesh[i], mesh[j]))
				return Validity::SelfIntersection;

	return Validity::OK;
}

TEST_CASE("IsValid") {
	dvec3 a(0, 0, 0);
	dvec3 b(1, 0, 0);
	dvec3 c(0, 1, 0);
	dvec3 d(0, 0, 1);
    REQUIRE(IsValid(std::vector<triangle3>{}) == Validity::TooFewFaces);
}

// Mesh properties
// ===============

// Volume of a valid polyhedron
real Volume(const Mesh3d& mesh) {
	real volume = 0;
	FOR_EACH(m, mesh) {
		real s = 0, z = 0;
		FOR_EACH_EDGE(a, b, m) {
			z += b->z;
			s += (a->y + b->y) * (a->x - b->x);
		}
		volume += z * s;
	}
	return volume / 6;
}

// TODO test with volume of cube (randomly rotated)

// Center of mass of a valid polyhedron
dvec3 CenterOfMass(const Mesh3d& mesh) {
	dvec3 P(0, 0, 0);
	real V = 0;
	FOR_EACH(f, mesh) {
		real v = dot(f.a, cross(f.b, f.c));
		P += (f.a + f.b + f.c) * v;
		V += v;
	}
	return P / (V * 4);
}

dmat3 full_mat3(real a) {
    const real m[] = { a, a, a, a, a, a, a, a, a};
    return glm::make_mat3(m);
}

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
dmat3 MomentOfInertia(const Mesh3d& mesh) {
	constexpr real a = 1 / 60., b = 1 / 120.;
    constexpr real canon[] = { a, b, b, b, a, b, b, b, a};
	const dmat3 canonical = glm::make_mat3(canon);
    dmat3 C = dmat3(0); // covariance
	FOR_EACH(f, mesh) {
		dmat3 A;
		FOR(i, 3)
			glm::row(A, i) = f[i]; // TODO: or column(i) = ?
		C += (transpose(A) * canonical * A) * determinant(A);
	}
	return full_mat3(C[0][0] + C[1][1] + C[2][2]) - C; // C -> I
}

// Is this valid mesh a convex polyhedron?
bool IsConvex(const Mesh3d& mesh) {
	// Extract all unique vertices
	std::unordered_set<dvec3> vertices;
	FOR_EACH(f, mesh)
		vertices.insert(f.begin(), f.begin() + 2);

	FOR_EACH(f, mesh) {
		// TODO can use the same trick from ConvexHull to avoid computing 1/sqrt()
		plane plane(f);
		FOR_EACH(v, vertices)
			if (plane.distance(v) > PlanarEpsilon)
				return false;
	}
	return true;
}

// Mesh generators
// ===============

// Returns empty vector if no solution (points are coplanar)
Mesh3d ConvexHullMesh(const std::vector<dvec3>& points) {
	// First two points (A and B) on the hull (extremes on X axis)
	size_t ai = 0, bi = 0;
	FOR(i, points.size()) {
		const dvec3& p = points[i];
		if (p.x < points[ai].x)
			ai = i;
		if (p.x > points[bi].x)
			bi = i;
	}
	if (ai == bi)
		return Mesh3d();
	dvec3 a = points[ai], b = points[bi];

	// Third point C on the hull (furthest from line AB)
	segment3 line(a, b);
	real max_dist2 = 0;
	size_t ci = 0;
	FOR(i, points.size()) {
		const dvec3& p = points[i];
		real dist2 = line_point_squared_distance(a, b, p);
		if (dist2 > max_dist2) {
			ci = i;
			max_dist2 = dist2;
		}
	}
	if (max_dist2 < squared(PlanarEpsilon))
		return Mesh3d();
	dvec3 c = points[ci];

	// Fourth point D on the hull (furthest from plane ABC)
	plane plane(a, b, c);
	real max_dist = 0;
	size_t di = 0;
	FOR(i, points.size()) {
		const dvec3& p = points[i];
		real dist = plane.distance(p);
		if (abs(dist) > abs(max_dist)) {
			di = i;
			max_dist = dist;
		}
	}
	if (abs(max_dist) < PlanarEpsilon)
		return Mesh3d();
	dvec3 d = points[di];

	// Construct initial tetrahedron hull,
	// All faces are oriented facing outside with right hand rule.
	Mesh3d hull;
	hull.reserve(4);
	if (max_dist < 0) {
		hull.push_back(triangle3(a, b, c));
		hull.push_back(triangle3(b, a, d));
		hull.push_back(triangle3(c, b, d));
		hull.push_back(triangle3(a, c, d));
	} else {
		hull.push_back(triangle3(c, b, a));
		hull.push_back(triangle3(a, b, d));
		hull.push_back(triangle3(b, c, d));
		hull.push_back(triangle3(c, a, d));
	}

	// Expand hull to include all remaining points outside of it
	FOR(pi, points.size()) {
		const dvec3& p = points[pi];
		if (pi == ai || pi == bi || pi == ci || pi == di)
			continue;

		// Remove faces on hull covered by new vertex
		std::unordered_set<segment3> open_edges;
		FOR(i, hull.size()) {
			const triangle3& f = hull[i];
			dvec3 normal = Normal(f.a, f.b, f.c);
			real dist = dot(normal, p - f.a);

			// If P is in front of face F
			if (dist > 0 && dist * dist > squared(PlanarEpsilon) * squared(normal)) {
				// Add edges of removed face to open_edges
				FOR_EACH_EDGE(a, b, f) {
					segment3 e(*a, *b);
					// If two faces that share an edge are both removed,
					// then their shared edge isn't open anymore.
					if (open_edges.erase(e.reverse()) == 0)
						open_edges.insert(e);
				}

				// Remove face
				hull[i] = hull[hull.size() - 1];
				hull.resize(hull.size() - 1);
				i -= 1;
			}
		}

		// For each open edge create a face that connects it with P
		FOR_EACH(e, open_edges)
			hull.push_back(triangle3(e.a, e.b, p));
	}
	return hull;
}

Mesh3d BoxMesh(real sx, real sy, real sz) {
	std::vector<dvec3> vertices;
	vertices.reserve(8);
	for (int x = -1; x <= 1; x += 2)
		for (int y = -1; y <= 1; y += 2)
			for (int z = -1; z <= 1; z += 2)
				vertices.push_back(dvec3(x * sx, y * sy, z * sz));
	return ConvexHullMesh(vertices);
}

Mesh3d SphereMesh(int vertices) {
	std::default_random_engine random;
	std::normal_distribution<real> gauss(0.0, 1.0);
	std::vector<dvec3> V(vertices);
	FOR(i, vertices)
		V[i] = glm::normalize(dvec3(gauss(random), gauss(random), gauss(random)));
	return ConvexHullMesh(V);
}

TEST_CASE("SphereMesh - ConvexHull - IsConvex") {
	Mesh3d mesh = SphereMesh(100);
	real v = 4.0 / 3 * M_PI;
	REQUIRE(abs(Volume(mesh) - v) / v < 0.15);
	REQUIRE(IsValid(mesh) == Validity::OK);
	REQUIRE(IsConvex(mesh));
}

void AddQuad(Mesh3d& mesh, dvec3 a, dvec3 b, dvec3 c, dvec3 d) {
	mesh.push_back(triangle3(a, b, c));
	mesh.push_back(triangle3(c, d, a));
}

Mesh3d CreateCrossMesh(real inner, real outer) {
	Mesh3d m;
	// TODO
	real a = inner, b = outer;

	// square face - Z axis
	AddQuad(m, dvec3(a, a, b), dvec3(-a, a, b), dvec3(-a, -a, b), dvec3(a, -a, b));
	AddQuad(m, dvec3(a, a, -b), dvec3(-a, -a, -b), dvec3(-a, a, -b), dvec3(a, -a, -b));
	// square face - Y axis
	AddQuad(m, dvec3(a, b, a), dvec3(-a, b, a), dvec3(-a, b, -a), dvec3(a, b, -a));
	AddQuad(m, dvec3(a, -b, a), dvec3(-a, -a, -a), dvec3(-a, -b, a), dvec3(a, -b, -a));
	// square face - X axis
	AddQuad(m, dvec3(b, a, a), dvec3(b, -a, a), dvec3(b, -a, -a), dvec3(b, a, -a));
	AddQuad(m, dvec3(-b, a, a), dvec3(-b, -a, -a), dvec3(-b, -a, a), dvec3(-b, a, -a));
	return m;
}

TEST_CASE("CreateCrossMesh") {
	Mesh3d cross = CreateCrossMesh(0.05, 0.25);
	//REQUIRE(IsValid(cross) == Validity::OK);
}

// Assumes center of sphere is at origin
real BoundingSphereRadius(const Mesh3d& mesh) {
	real squared_radius = 0;
	FOR_EACH(f, mesh)
		FOR(i, 3)
			squared_radius = std::max(squared_radius, squared(f[i]));
	return sqrt(squared_radius);
}

// Shape is rigid and immutable and purely geometric
class Shape {
public:
	Shape(const Mesh3d& mesh) {
		REQUIRE(IsValid(mesh) == Validity::OK); // TODO handle failure

		// Move to center of mass
		dvec3 center_of_mass = CenterOfMass(mesh);
		Mesh3d centered_mesh(mesh);
		FOR(i, mesh.size()) centered_mesh[i] -= center_of_mass;
		// TODO return back this center of mass, so objects isn't shifted

		m_volume = Volume(centered_mesh);
		// TODO m_inertia_tensor = MomentOfInertia(centered_mesh);
		// TODO rotate polyhedron, so inertia tensor only has primary axes (but return back its original orientation)
		// TODO so inertia tensor becomes inertia vector
		// TODO will this rotation make AABB smaller?

		m_mesh = std::move(centered_mesh);
		m_bounding_radius = BoundingSphereRadius(m_mesh);

		m_is_convex = ::IsConvex(m_mesh);
		// TODO init m_convex_planes

		// TODO init convex edges and vertices
		std::unordered_map<segment3, int> edge_to_face;
		FOR(i, m_mesh.size())
			FOR_EACH_EDGE(a, b, m_mesh[i])
				edge_to_face[segment3(*a, *b)] = i;

		std::unordered_set<dvec3> set_of_convex_vertices;
		FOR_EACH(face, m_mesh)
			FOR_EACH_EDGE(a, b, face)
				if (IsConcaveEdge(*a, *b)) {
					set_of_convex_vertices.erase(*a);
					set_of_convex_vertices.erase(*b);
				}
        FOR_EACH(v, set_of_convex_vertices)
		  m_convex_vertices.push_back(v);
	}

	const auto& Faces() const { return m_mesh; }
	const auto& ConvexEdges() const { return m_convex_edges; }
	const auto& ConvexVertices() const { return m_convex_vertices; }
	bool IsConvex() const { return m_is_convex; }
    bool IsConcaveEdge(dvec3, dvec3) const { return false; /*TODO*/ }

private:
	// TODO move faces and edges that collide more often to start of the vector

	Mesh3d m_mesh; // TODO to save memory: keep unique vertices in separate array and use 16bit pointers

	std::vector<segment3> m_convex_edges; // non-concave and non-planar
	std::vector<dvec3> m_convex_vertices; // non-concave and non-planar (needs to have more than a hemisphere around it open)

	dvec3 m_inertia_principal_axes; // assuming density of 1kg/m^3
	real m_volume;
	real m_bounding_radius;

	bool m_is_convex;
	std::vector<plane> m_convex_planes;

	// TODO is simple box mesh?
	// TODO maybe plane for each face? (or at least a unit normal)
	// TODO maybe bounding box: AABB or OBB?
	// TODO maybe reduced bounding convex hull that fully encloses mesh?
};

const auto random_directions = [](){
	std::array<dvec3, 128> dirs;
	std::default_random_engine rnd;
	std::normal_distribution<real> gauss(0.0, 1.0);
	FOR_EACH(d, dirs)
		d = glm::normalize(dvec3(gauss(rnd), gauss(rnd), gauss(rnd)));
	return dirs;
}();

real distance(dvec3 a, dvec3 b) { return l2Norm(a - b); }

real distance(dvec3 a, segment3 b) { return distance(a, b.nearest(a)); }
real distance(segment3 a, dvec3 b) { return distance(b, a); }

real distance(segment3 a, segment3 b) {
    return sqrt(segment3::squared_distance(a, b));
}

constexpr bool inside_triangle_prism(dvec3 p, triangle3 m, dvec3 normal) {
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + normal, p) > 0)
			return false;
	return true;
}

real distance(dvec3 p, triangle3 m) {
    dvec3 n = Normal(m);
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + n, p) > 0)
			return distance(p, segment3(*a, *b));
	return abs(dot(n, p - m.a));
}

real distance(triangle3 m, dvec3 v) { return distance(v, m); }

real distance(segment3 s, triangle3 m) {
    real d1 = segment3::squared_distance(s, segment3(m.a, m.b));
    real d2 = segment3::squared_distance(s, segment3(m.b, m.c));
    real d3 = segment3::squared_distance(s, segment3(m.c, m.a));
    real d = sqrt(min(d1, d2, d3));

    dvec3 n = Normal(m);
    if (inside_triangle_prism(s.a, m, n))
        d = std::min(d, abs(dot(n, s.a - m.a)));
    if (inside_triangle_prism(s.b, m, n))
        d = std::min(d, abs(dot(n, s.b - m.a)));
    return d;
}

real distance(triangle3 m, segment3 s) { return distance(s, m); }

real distance(triangle3 p, triangle3 q) {
    real d = std::numeric_limits<real>::max();
    // TODO only compute these 6 distances if inside triangle prisms
    FOR(i, 3)
        d = min(d, distance(p, q[i]), distance(p[i], q));
    // TODO optimize - accumulate with squared distances
	FOR_EACH_EDGE(pa, pb, p)
	   FOR_EACH_EDGE(qa, qb, q)
           d = std::min(d, distance(segment3(*pa, *pb), segment3(*qa, *qb)));
    return d;
}

struct transform3 {
	dmat3 orientation;
	dvec3 position;

    transform3(dvec3 position, dquat orientation)
        : orientation(orientation), position(position) { }

    dvec3 to_global(dvec3 v) const { return orientation * v + position; }
    dvec3 to_local(dvec3 v) const { return (v - position) * orientation; }

    transform3 inverse() {
        return *this; // TODO
    }
};

transform3 combine(transform3 a, transform3 b) {
    // TODO matrix multiply (with inverse)
    return a;
}

struct Body {
	Shape shape; // TODO make const

    // TODO: move these (and functions) into separate class
	dvec3 position;
	dquat orientation;
	dmat3 orientation_matrix;

	// from global frame to local frame
	dvec3 ToLocal(const dvec3& p) const {
		return (p - position) * orientation_matrix;
	}

	// from local frame to global frame
	dvec3 ToGlobal(const dvec3& p) const {
		return orientation_matrix * p + position;
	}

	segment3 ToLocal(const segment3& p) const {
		return segment3(ToLocal(p.a), ToLocal(p.b));
	}

	dvec3 ToGlobalDir(const dvec3& p) const {
		return orientation_matrix * p;
	}

	segment3 ToGlobal(const segment3& p) const {
		return segment3(ToGlobal(p.a), ToGlobal(p.b));
	}

	triangle3 ToLocal(const triangle3& p) const {
		return triangle3(ToLocal(p.a), ToLocal(p.b), ToLocal(p.c));
	}

	triangle3 ToGlobal(const triangle3& p) const {
		return triangle3(ToGlobal(p.a), ToGlobal(p.b), ToGlobal(p.c));
	}
};

bool IsVertexPenetratingBody(const dvec3& v_global, const Body& body) {
	// TODO bounding box / sphere check

	// TODO optimize for box
	// TODO optimize for convex mesh

	// TODO even if mesh is concave we can have a convex bounding shape to quickly test against
	// (only if convex hull has small number of faces)

	// TODO if possible to split concave mesh into small number of convex ones: do it
	//      (even splitting large concave into smaller concave helps)

	// TODO could be precomputed with a dense voxel table:
	// each voxel is one of: outside / convex maybe / concave maybe / inside
	// convex maybe has list of faces that intersect the voxel (1 face in best case)
	// in case of concave maybe run full ray intersect algorithm

	// TODO optimize convex planar faces with more than 3 vertices (ie. sides of a cube)

	dvec3 v = body.ToLocal(v_global);

	FOR_EACH(dir, random_directions) {
		int hits = 0;
		FOR_EACH(face, body.shape.Faces()) {
			// TODO TODO compare ray (v,dir) and face
			// if ray strictly goes through the triangle then hits ^= 1
			// else if ray intersects face plane in point outside of face then continue
			// else if ray is coplanar with face but far from face bounding circle then continue
			// else (ray hits edge, vertex or is coplanar with face) then restart while loop
		}
		return hits != 0;
	}
	cerr << "IsVertexPenetratingMesh failure" << endl;
	return false;
}

bool IsEdgePenetratingBody(const segment3& edge_global, const Body& body) {
	segment3 edge = body.ToLocal(edge_global);

	// TODO TODO At least one of these negative cases is needed for performance
	// TODO use voxel grid here faster positive and negative cases

	// TODO bounding box / sphere check
	// TODO optimize for box
	// TODO optimize for convex mesh

	real edge_length = l2Norm(edge.b - edge.a);
	// t_min and t_max needed to avoid edge endpoints being too close triangle an merely touching it
	real t_min = ContactEpsilon / edge_length;
	real t_max = 1 - t_min;

	FOR_EACH(face, body.shape.Faces()) {
		// compare edge and triangle
		dvec3 normal = Normal(face);
		real dd = dot(normal, edge.b - edge.a);
		// if edge is parallel to triangle plane
		constexpr real tiny = 1e-8;
		if (abs(dd) < tiny || abs(dd) < tiny * l2Norm(normal) * edge_length)
			continue;
		real t = dot(normal, face.a - edge.a) / dd;
		// if intersection of infinite edge with triangle plane is outside of edge interior
		if (t < t_min || t > t_max)
			continue;
		dvec3 w = edge.linear(t); // intersection of triangle plane and edge

		FOR_EACH_EDGE(a, b, face) {
			// if W is outside of triangle
			if (plane::sign(*a, *b, *a + normal, w) > 0)
				return false;
			// if edge is too close to triangle edge (could be contact instead)
			if (segment3::squared_distance(edge, segment3(*a, *b)) <= squared(ContactEpsilon))
				return false;
		}
		return true;
	}

	// Hard cases:
	// - edges don't go through interior of any triangles
	// - endpoints aren't inside polyhedron
	//   => endpoints can be outside and on the surface
	// segment can intersect surface at the vertex or interior of edge

	// Ignoring hard cases as they can't happen in simulator (only adversarial examples)

	// Idea: look at immediate surface around polyhedron vertex or edge that intersects segment
	// That surface can be convex
	return false;
}

struct OBB {
    dvec3 size;
};

pair<real, real> project_obb(dvec3 obb_position, dvec3 obb_size, dmat3 obb_orientation, dvec3 dir) {
    // TODO
	return pair<real, real>(0.0, 0.0);
}

bool intersects_obb(const Body& body_a, const Body& body_b) {
    // TODO
	return false;
}

bool AreSeparatedQuick(const Body& body_a, const Body& body_b) {
    // TODO sphere / sphere
    // TODO OBB / OBB
    // TODO convex / convex (if at least one of them is concave)
	return false;
}

bool AreBodiesPenetrating(const Body& body_a, const Body& body_b) {
	// Check all edges of one against the other
	FOR_EACH(edge_a, body_a.shape.ConvexEdges())
		if (IsEdgePenetratingBody(body_a.ToGlobal(edge_a), body_b))
				return true;
	FOR_EACH(edge_b, body_b.shape.ConvexEdges())
		if (IsEdgePenetratingBody(body_b.ToGlobal(edge_b), body_a))
				return true;
	// Also check a single vertex from one body against another, in case of
	// one body being completely inside the other.
	return IsVertexPenetratingBody(body_a.ToGlobal(body_a.shape.ConvexVertices()[0]), body_b)
		|| IsVertexPenetratingBody(body_b.ToGlobal(body_b.shape.ConvexVertices()[0]), body_a);
}

struct Contact {
	real squared_dist;
	dvec3 position; // mid point
	dvec3 normal; // unit normal from body B to body A
};

bool IsVertexTouchingTriangle(const dvec3& p, const triangle3& m, Contact& /*out*/contact) {
	dvec3 normal = Normal(m);

	// if P is outside of triangle
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + normal, p) > 0)
			return false;

	// P is inside of triangle
	dvec3 nearest = p - normal * dot(normal, p - m.a);
	contact.squared_dist = squared(nearest - p);
	if (contact.squared_dist > squared(ContactEpsilon))
		return false;
	contact.position = (nearest + p) / static_cast<real>(2);
	contact.normal = glm::normalize(normal);
	return true;
}

bool IsEdgeTouchingEdge(const segment3& p, const segment3& q, Contact& /*out*/contact) {
	dvec3 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	real aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr real tiny = 1e-8;

	real d = aa * bb - ab * ab;
	real s = ab * bc - bb * ac;
	real t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0)) {
		dvec3 P = p.a + A * (s / d);
		dvec3 Q = q.a + B * (t / d);
		contact.squared_dist = squared(P - Q);
		if (contact.squared_dist > squared(ContactEpsilon))
			return false;
		contact.position = (P + Q) / static_cast<real>(2);
		// Note: using AxB here instead of P-Q as P and Q will be very close.
		contact.normal = glm::normalize(cross(A, B));
		if (dot(P - Q, contact.normal) < 0)
			contact.normal = -contact.normal;
		return true;
	}

	return false;
}

// Assuming bodies aren't penetrating each other
std::vector<Contact> FindAllContacts(const Body& body_a, const Body& body_b) {
	std::vector<Contact> contacts;
	Contact contact;

	// vertex vs face contacts
	FOR_EACH(vertex_a, body_a.shape.ConvexVertices()) {
		dvec3 vertex_a_global = body_a.ToGlobal(vertex_a);
		// TODO bounding volume check?
		dvec3 vertex_a_local = body_b.ToLocal(vertex_a_global);
		FOR_EACH(face_b, body_b.shape.Faces())
			if (IsVertexTouchingTriangle(vertex_a_local, face_b, /*out*/contact)) {
				contact.position = body_b.ToGlobal(contact.position);
				contact.normal = body_b.ToGlobalDir(contact.normal);
				contacts.push_back(contact);
			}
	}
	FOR_EACH(vertex_b, body_b.shape.ConvexVertices()) {
		dvec3 vertex_b_global = body_b.ToGlobal(vertex_b);
		// TODO bounding volume check?
		dvec3 vertex_b_local = body_a.ToLocal(vertex_b_global);
		FOR_EACH(face_a, body_a.shape.Faces())
			if (IsVertexTouchingTriangle(vertex_b_local, face_a, /*out*/contact)) {
				contact.position = body_a.ToGlobal(contact.position);
				contact.normal = -body_a.ToGlobalDir(contact.normal);
				contacts.push_back(contact);
			}
	}

	// edge vs edge contacts
	FOR_EACH(edge_a, body_a.shape.ConvexEdges()) {
		segment3 edge_a_global = body_a.ToGlobal(edge_a);
		// TODO bounding volume check?
		FOR_EACH(edge_b, body_b.shape.ConvexEdges())
			if (IsEdgeTouchingEdge(edge_a_global, body_b.ToGlobal(edge_b), /*out*/contact))
				contacts.push_back(contact);
	}

	return contacts;
}

// returns 0 if interecting
real distance(const Body& body_a, const Body& body_b) {
    // TODO optimize
    // TODO need strict intersection (not just deep penetration)
    if (IsVertexPenetratingBody(body_a.ToGlobal(body_a.shape.ConvexVertices()[0]), body_b))
        return 0;

    real dist = 0;

	// vertex vs face
	FOR_EACH(vertex_a, body_a.shape.ConvexVertices()) {
		dvec3 vertex_a_local = body_b.ToLocal(body_a.ToGlobal(vertex_a));
		FOR_EACH(face_b, body_b.shape.Faces())
            dist = std::min(dist, distance(face_b, vertex_a_local));
	}
	FOR_EACH(vertex_b, body_b.shape.ConvexVertices()) {
		dvec3 vertex_b_local = body_a.ToLocal(body_b.ToGlobal(vertex_b));
		FOR_EACH(face_a, body_a.shape.Faces())
            dist = std::min(dist, distance(face_a, vertex_b_local));
	}

	// edge vs edge contacts
	FOR_EACH(edge_a, body_a.shape.ConvexEdges()) {
		segment3 edge_a_local = body_b.ToLocal(body_a.ToGlobal(edge_a));
		FOR_EACH(edge_b, body_b.shape.ConvexEdges())
            dist = std::min(dist, distance(edge_a_local, edge_b));
	}

    return dist;
}

// Dynamics and simulation
// =======================

struct DynamicBody : public Body {
	real mass;
	dvec3 velocity;
	dvec3 angular_velocity;
};

class Joint {
	// ball joint (common point) 3DF
	// hinge / axel joint (common edge) 1DF
	// cylindrical joint (common line) 2DF
	// prismatic joint 1DF (two common parallel lines)
};

class BallJoint {
	dvec3 pa, pb;
};

class HingeJoint {
	dvec3 pa, qa, pb, qb;
	dvec3 ra, rb; // just for computing relative angle (must be unit and normal to PQ)
};

// Open chain articulated only
// TODO describe children with relative coordinates, but allow computation of absolute ones
class ArticulatedBody {
	DynamicBody m_body;
	std::vector<pair<ArticulatedBody, Joint>> m_children;
};

template<typename T>
T EulersMethod(T x0, real t0, real dt, std::function<T(T x, real t)> derivative) {
	T k1 = dt * derivative(x0, t0);
	return x0 + k1;
}

template<typename T>
T MidpointMethod(T x0, real t0, real dt, std::function<T(T x, real t)> derivative) {
	T k1 = dt * derivative(x0, t0);
	T k2 = dt * derivative(x0 + k1 / 2, t0 + dt / 2);
	return x0 + k2;
}

// computes X[dt] given X[0] and derivative function of X at x0 and t
// T can be scalar, vector or matrix
template<typename T>
T RungeKutta4(T x0, real t0, real dt, std::function<T(T x, real d)> derivative) {
	T k1 = dt * derivative(x0, t0);
	T k2 = dt * derivative(x0 + k1 / 2, t0 + dt / 2);
	T k3 = dt * derivative(x0 + k2 / 2, t0 + dt / 2);
	T k4 = dt * derivative(x0 + k3, t0 + dt);
	return x0 + (k1 + 2 * (k2 + k3) + k4) / 6;
}

TEST_CASE("Numerical integrators") {
	// differential edquation: y' = y, y(0) = 1
	// exact solution: y(t) = e^t
	// simulating 320 steps from 0.0 to 4.0, step size 0.0125
	const auto& der = [](real e, real t) { return e; };
	real k = 4.0 / 320;

	real e = 1;
	FOR(i, 320) e = EulersMethod<real>(e, i * k, k, der);
	REQUIRE(abs(e - exp(4.0)) <= 1.34);

	e = 1;
	FOR(i, 320) e = MidpointMethod<real>(e, i * k, k, der);
	REQUIRE(abs(e - exp(4.0)) <= 0.0057);

	e = 1;
	FOR(i, 320) e = RungeKutta4<real>(e, i * k, k, der);
	REQUIRE(abs(e - exp(4.0)) <= 4.4e-08);
}
