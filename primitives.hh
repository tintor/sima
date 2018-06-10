#ifndef __PRIMITIVES_H__
#define __PRIMITIVES_H__

#include <cmath>
#include <random>

#include "tinyformat.h"
#include "for.hh"
#include "util.hh"

using real = double;

#define GLM_ENABLE_EXPERIMENTAL
#define GLM_FORCE_CXX14
#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtx/quaternion.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_access.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/gtx/transform.hpp"
using namespace glm;

static_assert(sizeof(dvec3) == sizeof(real) * 3);

inline bool lexicographical_less(dvec3 a, dvec3 b) {
	return a.x < b.x || (a.x == b.x && a.y < b.y) || (a.x == b.x && a.y == b.y && a.z < b.z);
}

namespace glm {
	inline std::ostream& operator<<(std::ostream& os, const dvec3& v) { return os << v.x << ' ' << v.y << ' ' << v.z; }
}

inline constexpr real squared(real a) { return a * a; }
inline real squared(dvec3 a) { return dot(a, a); }
inline real squared(vec4 a) { return dot(a, a); }

inline constexpr long double operator "" _deg(long double a) { return a * (M_PI / 180); }
inline constexpr long double operator "" _deg(unsigned long long a) { return a * (M_PI / 180); }

inline std::string deg(long double r) { return tfm::format("%llg deg", r * (180 / M_PI)); }

inline auto angle(dvec3 a, dvec3 b) {
	// atan2 is numericaly better than acos when angle is very small

	// make the vectors equal length first (without division)
	real aa = l2Norm(a), bb = l2Norm(b);
	a *= bb;
	b *= aa;
	return 2 * atan2(l2Norm(a - b), l2Norm(a + b));
}

inline constexpr real clamp(real t, real min = 0, real max = 1) {
	if (t < min) return min;
	if (t > max) return max;
	return t;
}

// infinite on both sides unlike segment3
struct line3 {
	dvec3 a, b;
};

struct ray3
{
    // TODO is dir constrained to UnitVector?
    dvec3 origin, dir;
};

// infinite line vs. point
real line_point_squared_distance(dvec3 a, dvec3 b, dvec3 p);

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

	// Return point on segment that is nearest to P
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
	static std::pair<segment3, NearestCase> Nearest(segment3 p, segment3 q);

	static real squared_distance(segment3 p, segment3 q);
};

inline auto angle(segment3 p, segment3 q) {
	return angle(p.b - p.a, q.b - q.a);
}

// uniform inside a cube
template<typename RandomEngine>
dvec3 random_vector(RandomEngine& rnd, real a = -1.0, real b = 1.0) {
	std::uniform_real_distribution<real> dist(a, b);
	return dvec3(dist(rnd), dist(rnd), dist(rnd));
}

// uniform on the unit sphere surface
template<typename RandomEngine>
dvec3 random_direction(RandomEngine& rnd) {
	std::normal_distribution<real> dist(0.0, 1.0);
	return normalize(dvec3(dist(rnd), dist(rnd), dist(rnd)));
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

inline dvec3 Normal(const dvec3& a, const dvec3& b, const dvec3& c) {
	return glm::cross(b - a, c - a);
}

inline dvec3 Normal(const triangle3& p) {
	return Normal(p.a, p.b, p.c);
}

// How close a point needs to be to a plane to be considered on the plane?
constexpr real PlanarEpsilon = 1e-6;

// How close two features need to be to be considered touching?
// TODO can it be same as PlanarEpsilon?
constexpr real ContactEpsilon = 10 * PlanarEpsilon;

struct plane {
	dvec3 normal; // must be unit vector
	real d;

	plane() { }

	plane(const dvec3& a, const dvec3& b, const dvec3& c)
		: normal(normalize(Normal(a, b, c))), d(dot(-normal, a)) { }

	plane(const triangle3& m) : plane(m.a, m.b, m.c) {}

	// Signed distance!
	real distance(const dvec3& p) const {
		return dot(normal, p) + d;
	}

	// > 0, if P is on the positive side of plane ABC (right hand rule)
	static real sign(const dvec3& a, const dvec3& b, const dvec3& c, const dvec3& p) {
		// TODO this is just a mixed product - which is 3x3 determinant?
		return dot(cross(b - a, c - a), p - a);
	}
};

inline constexpr bool intersects(const segment3& q, const plane& p) {
	real a = p.distance(q.a);
	real b = p.distance(q.b);
	real e = PlanarEpsilon;
	return (a <= e || b <= e) && (a >= -e || b >= -e);
}

inline constexpr bool intersects(const triangle3& q, const plane& p) {
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

inline dvec3 mini(const triangle3& p) {
	dvec3 e;
	FOR(i, 3)
		e[i] = min(p.a[i], p.b[i], p.c[i]);
	return e;
}

inline dvec3 maxi(const triangle3& p) {
	dvec3 e;
	FOR(i, 3)
		e[i] = max(p.a[i], p.b[i], p.c[i]);
	return e;
}

inline bool DisjointIntervals(real amin, real amax, real bmin, real bmax) {
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

template<typename T>
void hash_combine(size_t& seed, const T& value) {
	// The code is from `hash_combine` function of the Boost library. See
	// http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine .
	seed ^= std::hash<T>()(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
	template<> struct hash<dvec3> {
		size_t operator()(const dvec3& a) const {
			size_t seed = 0;
			FOR(i, 3)
				hash_combine(seed, a[i]);
			return seed;
		}
    };

  template <typename A, typename B> struct hash<pair<A, B>> {
    size_t operator()(const pair<A, B>& x) const {
		size_t seed = 0;
		hash_combine(seed, x.first);
		hash_combine(seed, x.second);
	    return seed;
    }
  };

  template<> struct hash<segment3> {
    size_t operator()(const segment3& x) const {
		size_t seed = 0;
		hash_combine(seed, x.a);
		hash_combine(seed, x.b);
	    return seed;
    }
  };
}

inline dmat3 full_mat3(real a) {
    const real m[] = { a, a, a, a, a, a, a, a, a};
    return glm::make_mat3(m);
}

real distance(const dvec3& a, const dvec3& b);
real distance(const dvec3& a, const segment3& b);
real distance(const segment3& a, const dvec3& b);
real distance(const segment3& a, const segment3& b);

constexpr bool inside_triangle_prism(const dvec3& p, const triangle3& m, const dvec3& normal);

real distance(const dvec3& p, const triangle3& m);
real distance(const dvec3& v, const triangle3& m, const plane& p);
real distance(const triangle3& m, const dvec3& v);

bool intersects(const line3& e, const triangle3& m);
bool intersects2(const line3& e, const triangle3& m);
bool intersection(const line3& e, const triangle3& m, /*out*/vec3& result);
bool intersects_in_point(const segment3& e, const triangle3& m);
bool intersects_in_point(const ray3& e, const triangle3& m);

real disjoint_distance(const segment3& e, const triangle3& m);
real distance(const segment3& e, const triangle3& m);

real distance(const triangle3 m, const segment3& e);
real disjoint_distance(const triangle3& p, const triangle3& q);
real distance(const triangle3& p, const triangle3& q);

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

transform3 combine(const transform3& a, const transform3& b);

// Angle between oriented triangles ABC and BAD.
// If edge is convex then angle will be <PI
// If edge is planar then angle will be =PI
// If edge is concave than angle will be >PI
real edge_angle(const dvec3& a, const dvec3& b, const dvec3& c, const dvec3& d);

struct sphere {
	dvec3 center;
	real radius;
};

sphere merge_spheres(const sphere& a, const sphere& b);

#endif
