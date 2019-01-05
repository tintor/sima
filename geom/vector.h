#pragma once
#include <core/format.h>
#include <geom/simd.h>
#include <core/util.h>
#include <core/align_alloc.h>
#include <random>

namespace std {

template<>
class allocator<double3> : public AlignAlloc<double3> { };
template<>
class allocator<double4> : public AlignAlloc<double4> { };

}

#define REQUIRE_NEAR(A, B, T) do { auto aa = A; auto bb = B; auto tt = T; auto dd = length(aa - bb); ASSERT_ALWAYS(dd <= tt, "a = %s\nb = %s\nexpected <= %s, actual = %s", aa, bb, tt, dd); } while (false);

#define REQUIRE_EQUAL(A, B) do { auto aa = A; auto bb = B; ASSERT_ALWAYS(all(aa == bb), "a = %s\nb = %s\n|a-b| = %g", aa, bb, abs(aa - bb)); } while (false);

inline double2 d2(double x, double y) { return {x, y}; }
inline double3 d3(double x, double y, double z) { return {x, y, z}; }
inline double4 d4(double x, double y, double z, double w) { return {x, y, z, w}; }

inline double3 extend(double2 v, double z) { return {v.x, v.y, z}; }
inline double4 extend(double3 v, double w) { return {v.x, v.y, v.z, w}; }

inline long2 lbroad2(long a) { return {a, a}; }
inline long4 lbroad4(long a) { return {a, a, a, a}; }

inline double2 broad2(double a) { return {a, a}; }
inline double3 broad3(double a) { return {a, a, a}; }
inline double4 broad4(double a) { return {a, a, a, a}; }

inline void broadcast(double2& a, double b) { a = {b, b}; }
inline void broadcast(double3& a, double b) { a = {b, b, b}; }
inline void broadcast(double4& a, double b) { a = {b, b, b, b}; }

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

inline bool lex_less(double4 a, double4 b) {
	if (a.x < b.x) return true;
	if (a.x > b.x) return false;
	if (a.y < b.y) return true;
	if (a.y > b.y) return false;
	if (a.z < b.z) return true;
	if (a.z > b.z) return false;
	return a.w < b.w;
}

inline double2 sqrt(double2 a) { return _mm_sqrt_pd(a); }
inline double3 sqrt(double3 a) { return double4(_mm256_sqrt_pd(cast4(a))).xyz; }
inline double4 sqrt(double4 a) { return _mm256_sqrt_pd(a); }

// returns 1.0 or -1.0 for each component (depending if >= +0, or <= -0)
inline double4 sign_no_zero(double4 d) {
	return bit_and(bit_or(d, broad4(1.0)), broad4(-1.0));
}
inline double3 sign_no_zero(double3 d) {
	return sign_no_zero(cast4(d)).xyz;
}

inline float4 vdot(float4 a, float4 b) { return _mm_dp_ps(a, b, 0b11111111); }

inline double2 vdot(double2 a, double2 b) { return _mm_dp_pd(a, b, 0b11111111); }
// 4 instructions vs 1 instruction for float4
inline double4 vdot(double4 a, double4 b) {
	// There is no _mm256_dp_pd! But there is _mm_dp_ps for floats!
	// return _mm256_dp_pd(a, a, 0b11111111);
	double4 ab = a * b;
	double4 q = _mm256_hadd_pd(ab, ab);
	return q + q.zwxy;
}

inline float dot(float4 a, float4 b) { return vdot(a, b).x; }

inline double dot(double2 a, double2 b) { return vdot(a, b).x; }
inline double dot(double3 a, double3 b) { double3 q = a * b; return q.x + q.y + q.z; }
inline double dot(double4 a, double4 b) { return vdot(a, b).x; }

inline constexpr double squared(double a) { return a * a; }
inline double squared(double2 a) { return dot(a, a); }
inline double squared(double3 a) { return dot(a, a); }
inline double squared(double4 a) { return dot(a, a); }

inline double2 vlength(double2 a) { return sqrt(vdot(a, a)); }
inline double4 vlength(double4 a) { return sqrt(vdot(a, a)); }

inline constexpr double length(double a) { return (a >= 0) ? a : -a; }
inline double length(double2 a) { return vlength(a).x; }
inline double length(double3 a) { return sqrt(dot(a, a)); }
inline double length(double4 a) { return vlength(a).x; }

// TODO div(sqrt) might be slower than mul(rsqrt)
inline auto normalize(double2 a) { return a / vlength(a); }
inline auto normalize(double3 a) { return a / length(a); }
inline auto normalize(double4 a) { return a / vlength(a); }

template<typename V> bool is_unit(V v) { return ordered(1-1e-12, squared(v), 1+1e-12); }

template<int M> double2 round(double2 v, int mode) { return _mm_round_pd(v, M | _MM_FROUND_NO_EXC); }
template<int M> double3 round(double3 v, int mode) { double4 e = _mm256_round_pd(cast4(v), M | _MM_FROUND_NO_EXC); return e.xyz; }
template<int M> double4 round(double4 v, int mode) { return _mm256_round_pd(v, M | _MM_FROUND_NO_EXC); }

template<typename V> V round_nearest(V v) { return round<_MM_FROUND_TO_NEAREST_INT>(v); }
template<typename V> V round_zero(V v) { return round<_MM_FROUND_TO_ZERO>(v); }
template<typename V> V round_up(V v) { return round<_MM_FROUND_TO_POS_INF>(v); }
template<typename V> V round_down(V v) { return round<_MM_FROUND_TO_NEG_INF>(v); }

constexpr ulong SignBit64 = ulong(1) << 63;
inline double2 abs(double2 v) { return bit_and(v, lbroad2(~SignBit64)); }
inline double3 abs(double3 v) { double4 e = bit_and(cast4(v), lbroad4(~SignBit64)); return e.xyz; }
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

inline double3 cross(double3 a, double3 b) {
	return {a.y * b.z - b.y * a.z,
            a.z * b.x - b.z * a.x,
            a.x * b.y - b.x * a.y};
}

// returns angle in range [0, PI]
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

// det(A) == det(transposed(A))
inline double det(double2 a, double2 b) {
	return a.x * b.y - b.x * a.y;
}

// det(A) == det(transposed(A))
// 3x3 det
inline double det(double3 a, double3 b, double3 c) {
	return a.x * det(b.yz, c.yz)
         - b.x * det(a.yz, c.yz)
		 + c.x * det(a.yz, b.yz);
}

// det(A) == det(transposed(A))
inline double det(double4 a, double4 b, double4 c, double4 d) {
	return a.x * det(b.yzw, c.yzw, d.yzw)
         - b.x * det(a.xzw, c.xzw, d.xzw)
		 + c.x * det(a.xyw, b.xyw, d.xyw)
		 - d.x * det(a.xyz, b.xyz, c.xyz);
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

// uniform point inside a cube
template<typename RND>
double3 uniform3(RND& rnd, double min, double max) {
	return {uniform(rnd, min, max), uniform(rnd, min, max), uniform(rnd, min, max)};
}

template<typename RND>
double4 uniform4(RND& rnd, double min, double max) {
	return {uniform(rnd, min, max), uniform(rnd, min, max), uniform(rnd, min, max), uniform(rnd, min, max)};
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

template<typename RND>
double4 normal4(RND& rnd, double mean, double stdev) {
	return {normal(rnd, mean, stdev), normal(rnd, mean, stdev), normal(rnd, mean, stdev), normal(rnd, mean, stdev)};
}

template<typename RND>
double2 uniform_dir2(RND& rnd) {
	double a = uniform(rnd, 0, 2 * PI);
	return {cos(a), sin(a)};
}

template<typename RND>
double3 uniform_dir3(RND& rnd) { return normalize(normal3(rnd, 0, 1)); }

template<typename RND>
double4 uniform_dir4(RND& rnd) { return normalize(normal4(rnd, 0, 1)); }

inline bool colinear(double3 a, double3 b, double3 c) {
	return squared(cross(b - a, c - a)) <= squared(1e-6);
}

inline double3 avg(double3 a, double3 b) { return (a + b) / 2; }

inline double3 compute_normal(double3 a, double3 b, double3 c) { return cross(b - a, c - a); }

inline void remove_dups(vector<double3>& vertices) {
	::remove_dups(vertices,
		static_cast<bool(*)(double3, double3)>(&lex_less),
		static_cast<bool(*)(double3, double3)>(&equal));
}
