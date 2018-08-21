#pragma once
#include "format.h"
#include "vec.h"
#include "util.h"
#include <random>

inline long2 lbroad2(long a) { return {a, a}; }
inline double2 broad2(double a) { return {a, a}; }
inline long4 lbroad4(long a) { return {a, a, a, a}; }
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

inline double2 sqrt(double2 a) { return _mm_sqrt_pd(a); }
inline double4 sqrt(double4 a) { return _mm256_sqrt_pd(a); }
inline __m128d bit_and(__m128d a, __m128d b) { return _mm_and_pd(a, b); }
inline __m256d bit_and(__m256d a, __m256d b) { return _mm256_and_pd(a, b); }

inline float4 vdot(float4 a, float4 b) { return _mm_dp_ps(a, b, 0b11111111); }

inline double2 vdot(double2 a, double2 b) { return _mm_dp_pd(a, b, 0b11111111); }
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

inline float dot(float4 a, float4 b) { return vdot(a, b).x; }
inline double dot(double2 a, double2 b) { return vdot(a, b).x; }
inline double dot(double3 a, double3 b) { return vdot(a, b).x; }
inline double dot(double4 a, double4 b) { return vdot(a, b).x; }

inline constexpr auto squared(double a) { return a * a; }

inline double squared(double2 a) { return dot(a, a); }
inline double squared(double3 a) { return dot(a, a); }
inline double squared(double4 a) { return dot(a, a); }

inline double2 vlength(double2 a) { return sqrt(vdot(a, a)); }
inline double4 vlength(double3 a) { return sqrt(vdot(a, a)); }
inline double4 vlength(double4 a) { return sqrt(vdot(a, a)); }
template<typename V> double length(V a) { return vlength(a).x; }

// TODO div(sqrt) might be slower than mul(rsqrt)
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

inline double3 cross(double3 a, double3 b) {
	return {a.y * b.z - b.y * a.z,
            a.z * b.x - b.z * a.x,
            a.x * b.y - b.x * a.y};
}

inline double4 cross(double4 a, double4 b) {
	assert(a.w == 0);
	assert(b.w == 0);
	return {a.y * b.z - b.y * a.z,
            a.z * b.x - b.z * a.x,
            a.x * b.y - b.x * a.y,
			0};
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

// det(A) == det(transposed(A))
inline double det(double2 a, double2 b) {
	return a.x * b.y - b.x * a.y;
}

// det(A) == det(transposed(A))
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

// TODO define struct vec3
// w can always be assumed to be 0
// with alignas(32)
// operator ==
// hash

struct alignas(32) vec3 {
private:
	double4 s;
};

struct alignas(32) point3 {
	point3() : s{NAN, NAN, NAN, 1} { }
	point3(double3 v) : s{v.x, v.y, v.z, 1} { }
	point3(double x, double y, double z) : s{x, y, z, 1} { }
	point3(const point3& p) : s(p.s) { }

	bool operator==(point3 v) const { return equal(s, v.s); }
	bool operator!=(point3 v) const { return !equal(s, v.s); }

	bool operator<(point3 v) const {
		if (s.x < v.s.x) return true;
		if (s.x > v.s.x) return false;
		if (s.y < v.s.y) return true;
		if (s.y > v.s.y) return false;
		return s.z < v.s.z;
	}

	double4 operator-(point3 v) const { return s - v.s; }

	point3 operator+(double4 v) const { return s + v; }
	point3 operator-(double4 v) const { return s - v; }
	point3 operator+(double3 v) const { return s + double4{v.x, v.y, v.z, 0}; }
	point3 operator-(double3 v) const { return s - double4{v.x, v.y, v.z, 0}; }

	point3& operator+=(double4 v) { s += v; return *this; }
	point3& operator-=(double4 v) { s -= v; return *this; }
	point3& operator+=(double3 v) { s += {v.x, v.y, v.z, 0}; return *this; }
	point3& operator-=(double3 v) { s -= {v.x, v.y, v.z, 0}; return *this; }

	double x() const { return s.x; }
	double y() const { return s.y; }
	double z() const { return s.z; }

	friend void broadcast(point3& a, double b);
	friend point3 avg(point3 a, point3 b);
	friend point3 vmin(point3 a, point3 b);
	friend point3 vmax(point3 a, point3 b);
	friend double dot(point3 a, double4 b);
	friend double dot(double4 a, point3 b);
	friend float dot(point3 a, float4 b);
	friend float dot(float4 a, point3 b);

	operator double4() { return s; }
private:
	point3(double4 v) : s(v) { }
	double4 s;
};

inline void format_e(string& s, string_view spec, point3 v) {
	format_e(s, spec, v.x());
	s += ' ';
	format_e(s, spec, v.y());
	s += ' ';
	format_s(s, spec, v.z());
}

namespace std {

template<> struct hash<point3> {
	size_t operator()(point3 v) const { return ::hash(v.x(), v.y(), v.z()); }
};

}

inline void broadcast(point3& a, double b) { broadcast(a.s, b); }
inline point3 avg(point3 a, point3 b) { return (a.s + b.s) / 2; }
inline point3 vmin(point3 a, point3 b) { return vmin(a.s, b.s); }
inline point3 vmax(point3 a, point3 b) { return vmax(a.s, b.s); }
inline double dot(point3 a, double4 b) { return dot(a.s, b); }
inline double dot(double4 a, point3 b) { return dot(a, b.s); }
inline float dot(point3 a, float4 b) { return dot(vconvert(a.s, float4), b); }
inline float dot(float4 a, point3 b) { return dot(a, vconvert(b.s, float4)); }

inline double3 compute_normal(double3 a, double3 b, double3 c) { return cross(b - a, c - a); }
inline double4 compute_normal(point3 a, point3 b, point3 c) { return cross(b - a, c - a); }
