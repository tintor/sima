#pragma once
#include "format.h"
#include "vec.h"
#include "util.h"
#include <random>

// TODO different classes for point2/point3 and vec2/vec3
// point-point -> vec
// point+point error!
// point*k error!
// point+-vec -> point
// vec*k -> vec
// vec+-vec -> vec
// explicit convert point <-> vec
// point3 is double4 with w=1
// vec3 is double4 with w=0

// TODO maybe just point class with w=1 and operator==?
struct point3 {
	double4 s;

	point3() : s{NAN, NAN, NAN, 1} { }
	point3(double x, double y, double z) : s{x, y, z, 1} { }
	point3(const point3& p) : s(p.s) { }
};

/*struct vec4 {
	double4 s;

	vec4() : s{NAN, NAN, NAN, NAN} { }
	vec4(double x, double y, double z, double w=0) : s{x, y, z, w} { }
	vec4(const vec4& v) : s(v.s) { }
};*/

// TODO define proper vec2 and vec3 on top of double2 and double3
// - auto conversions
// - natural == operator
// - can use both {} and () initializer

// SIMD helpers
inline double3 d3(double4 v) { return v.xyz; }
inline double4 d4(double3 v) { return vshuffle(v, v, 0, 1, 2, -1); }
inline long2 lbroad2(long a) { return {a, a}; }
inline double2 broad2(double a) { return {a, a}; }
inline long4 lbroad4(long a) { return {a, a, a, a}; }
inline double4 broad4(double a) { return {a, a, a, a}; }

inline bool lex_less(double2 a, double2 b) {
	if (a.x < b.x) return true;
	if (a.x > b.x) return false;
	return a.y < b.y;
}

inline bool lex_less(double3 a, double3 b) {
	if (a.x < b.x) return true;
	if (a.x > b.x) return false;
	if (a.y < b.y) return true;
	if (a.y > b.y) return false;
	return a.z < b.z;
}

inline double2 sqrt(double2 a) { return _mm_sqrt_pd(a); }
inline double4 sqrt(double4 a) { return _mm256_sqrt_pd(a); }
inline __m128d bit_and(__m128d a, __m128d b) { return _mm_and_pd(a, b); }
inline __m256d bit_and(__m256d a, __m256d b) { return _mm256_and_pd(a, b); }

inline double2 vdot(double2 a, double2 b) { return _mm_dp_pd(a, a, 0b11111111); }
// 4 instructions vs 1 instruction for float3
inline double4 vdot(double3 a, double3 b) {
	// There is no _mm256_dp_pd!
	// return _mm256_dp_pd(d4(a), d4(a), 0b01111111);
	double4 ab = d4(a) * d4(b);
	assert(ab.w == 0); // convention that w=0 for all 3d points
	double4 q = _mm256_hadd_pd(ab, ab);
	return q + q.zwxy;
}
// 4 instructions vs 1 instruction for float4
inline double4 vdot(double4 a, double4 b) {
	// There is no _mm256_dp_pd! But there is _mm_dp_ps for floats!
	// return _mm256_dp_pd(a, a, 0b11111111);
	double4 ab = a * b;
	double4 q = _mm256_hadd_pd(ab, ab);
	return q + q.zwxy;
}

template<typename V> auto dot(V a, V b) { return vdot(a, b).x; }

inline constexpr auto squared(double a) { return a * a; }
template<typename V> auto squared(V a) { return dot(a, a); }

inline double2 vlength(double2 a) { return sqrt(vdot(a, a)); }
inline double4 vlength(double3 a) { return sqrt(vdot(a, a)); }
inline double4 vlength(double4 a) { return sqrt(vdot(a, a)); }
template<typename V> double length(V a) { return vlength(a).x; }

inline auto normalize(double2 a) { return a / vlength(a); }
inline auto normalize(double3 a) { return d3(d4(a) / vlength(a)); }
inline auto normalize(double4 a) { return a / vlength(a); }

template<typename V> bool is_unit(V v) { return ordered(1-1e-12, squared(v), 1+1e-12); }

template<int M> double2 round(double2 v, int mode) { return _mm_round_pd(v, M | _MM_FROUND_NO_EXC); }
template<int M> double3 round(double3 v, int mode) { return d3(_mm256_round_pd(d4(v), M | _MM_FROUND_NO_EXC)); }
template<int M> double4 round(double4 v, int mode) { return _mm256_round_pd(v, M | _MM_FROUND_NO_EXC); }

template<typename V> V round_nearest(V v) { return round<_MM_FROUND_TO_NEAREST_INT>(v); }
template<typename V> V round_zero(V v) { return round<_MM_FROUND_TO_ZERO>(v); }
template<typename V> V round_up(V v) { return round<_MM_FROUND_TO_POS_INF>(v); }
template<typename V> V round_down(V v) { return round<_MM_FROUND_TO_NEG_INF>(v); }

constexpr ulong SignBit64 = ulong(1) << 63;
inline double2 abs(double2 v) { return bit_and(v, lbroad2(~SignBit64)); }
inline double3 abs(double3 v) { return d3(bit_and(d4(v), lbroad4(~SignBit64))); }
inline double4 abs(double4 v) { return bit_and(v, lbroad4(~SignBit64)); }

inline double3 any_normal(double3 v) {
	double3 a = abs(v);
	if (a.x <= a.y && a.x <= a.z)
		return {0, -v.z, v.y};
	if (a.y <= a.z)
		return {-v.z, 0, v.x};
	return {-v.y, v.x, 0};
}

inline double cross(double2 a, double2 b) { return a.x * b.y - b.x * a.y; }

inline double3 cross(double3 a, double3 b){
	return {a.y * b.z - b.y * a.z,
            a.z * b.x - b.z * a.x,
            a.x * b.y - b.x * a.y};
}

inline double3 compute_normal(double3 a, double3 b, double3 c) {
	return cross(b - a, c - a);
}

// returns angle in range [0, PI)
inline auto angle(double3 a, double3 b) {
	return std::atan2(length(cross(a, b)), dot(a, b));
}

// solid angle between triangle and origin
inline double solid_angle(double3 A, double3 B, double3 C) {
    double y = dot(A, cross(B, C));
    double a = length(A), b = length(B), c = length(C);
    double x = a * b * c + c * dot(A, B) + b * dot(A, C) + a * dot(B, C);
    return 2 * std::atan2(y, x);
}

template<typename RND>
double uniform(RND& rnd, double min, double max) {
	return std::uniform_real_distribution<double>(min, max)(rnd);
}

// uniform inside a cube
template<typename RND>
double2 uniform2(RND& rnd, double min, double max) {
	return {uniform(rnd, min, max), uniform(rnd, min, max)};
}

// uniform inside a cube
template<typename RND>
double3 uniform3(RND& rnd, double min, double max) {
	return {uniform(rnd, min, max), uniform(rnd, min, max), uniform(rnd, min, max)};
}

template<typename RND>
inline double normal(RND& rnd, double mean, double stdev) {
	return std::normal_distribution<double>(mean, stdev)(rnd);
}

template<typename RND>
double2 normal2(RND& rnd, double mean, double stdev) {
	return {normal(rnd, mean, stdev), normal(rnd, mean, stdev)};
}

template<typename RND>
double3 normal3(RND& rnd, double mean, double stdev) {
	return {normal(rnd, mean, stdev), normal(rnd, mean, stdev), normal(rnd, mean, stdev)};
}

// uniform on the unit sphere surface
template<typename RND>
double3 uniform_dir3(RND& rnd) { return normalize(normal3(rnd, 0, 1)); }

inline bool colinear(double3 a, double3 b, double3 c) {
	return squared(cross(b - a, c - a)) <= squared(1e-6);
}
