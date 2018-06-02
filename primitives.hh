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

	dvec3 nearest_to_line(dvec3 p) const {
		dvec3 d = b - a;
		return a + d * (dot(p - a, d) / dot(d, d));
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

// 72 bytes (for doubles)
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

	real squared_area_x4() const {
		return squared(cross(b - a, c - a));
	}

	real area() const {
		return sqrt(squared_area_x4()) / 2;
	}
};

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
	uint32_t m_rank;
};

namespace std {
	// Hash function for matrix and vector.
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

dmat3 full_mat3(real a) {
    const real m[] = { a, a, a, a, a, a, a, a, a};
    return glm::make_mat3(m);
}

template<typename RND>
dvec3 random_direction(RND& rnd) {
	std::normal_distribution<real> gauss(0.0, 1.0);
	return glm::normalize(dvec3(gauss(rnd), gauss(rnd), gauss(rnd)));
}

real distance(const dvec3& a, const dvec3& b) { return l2Norm(a - b); }

real distance(const dvec3& a, const segment3& b) { return distance(a, b.nearest(a)); }
real distance(const segment3& a, const dvec3& b) { return distance(b, a); }

real distance(const segment3& a, const segment3& b) {
    return sqrt(segment3::squared_distance(a, b));
}

constexpr bool inside_triangle_prism(const dvec3& p, const triangle3& m, const dvec3& normal) {
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + normal, p) > 0)
			return false;
	return true;
}

real distance(const dvec3& p, const triangle3& m) {
    dvec3 n = Normal(m);
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + n, p) > 0)
			return distance(p, segment3(*a, *b));
	return abs(dot(n, p - m.a));
}

real distance(const triangle3& m, const dvec3& v) { return distance(v, m); }

real distance(const segment3& s, const triangle3& m) {
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

real distance(const triangle3 m, const segment3& s) { return distance(s, m); }

real distance(const triangle3& p, const triangle3& q) {
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

struct transform3 {
	dmat3 orientation;
	dvec3 position;

    transform3(const dvec3& position, dquat orientation)
        : orientation(orientation), position(position) { }

    dvec3 to_global(dvec3 v) const { return orientation * v + position; }
    dvec3 to_local(dvec3 v) const { return (v - position) * orientation; }

    segment3 to_local(const segment3& p) const {
		return segment3(to_local(p.a), to_local(p.b));
	}

	dvec3 to_global_dir(const dvec3& p) const {
		return orientation * p;
	}

	segment3 to_global(const segment3& p) const {
		return segment3(to_global(p.a), to_global(p.b));
	}

	triangle3 to_local(const triangle3& p) const {
		return triangle3(to_local(p.a), to_local(p.b), to_local(p.c));
	}

	triangle3 to_global(const triangle3& p) const {
		return triangle3(to_global(p.a), to_global(p.b), to_global(p.c));
	}
};

transform3 combine(const transform3& a, const transform3& b) {
    // TODO create two dmat4 from two transform3
	dmat4 a4, b4;
    // TODO figure out combination formula: ie: a * inverse(b)
	dmat4 c4 = a4 * inverse(b4);
    // TODO convert c4 back to transform3

    return a;
}

// Angle between oriented triangles ABC and BAD.
// If edge is convex then angle will be <PI
// If edge is planar then angle will be =PI
// If edge is concave than angle will be >PI
real edge_angle(const dvec3& a, const dvec3& b, const dvec3& c, const dvec3& d) {
	// TODO
	return 0;
}
