#pragma once
#include "format.h"
#include "vec.h"
#include "util.h"
#include "align_alloc.h"
#include <random>

inline long2 lbroad2(long a) { return {a, a}; }
inline double2 broad2(double a) { return {a, a}; }
inline long4 lbroad4(long a) { return {a, a, a, a}; }
inline double4 broad4(double a) { return {a, a, a, a}; }

inline void broadcast(double2& a, double b) { a = {b, b}; }
inline void broadcast(double4& a, double b) { a = {b, b, b, b}; }

inline bool lex_less(double2 a, double2 b) {
	if (a.x < b.x) return true;
	if (a.x > b.x) return false;
	return a.y < b.y;
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
inline double4 sqrt(double4 a) { return _mm256_sqrt_pd(a); }
inline __m128d bit_and(__m128d a, __m128d b) { return _mm_and_pd(a, b); }
inline __m256d bit_and(__m256d a, __m256d b) { return _mm256_and_pd(a, b); }
inline __m128d bit_or(__m128d a, __m128d b) { return _mm_or_pd(a, b); }
inline __m256d bit_or(__m256d a, __m256d b) { return _mm256_or_pd(a, b); }

// returns 1.0 or -1.0 for each component (depending if >= +0, or <= -0)
inline double4 sign_no_zero(double4 d) {
	return bit_and(bit_or(d, broad4(1.0)), broad4(-1.0));
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
inline double dot(double4 a, double4 b) { return vdot(a, b).x; }

inline constexpr auto squared(double a) { return a * a; }

inline double squared(double2 a) { return dot(a, a); }
inline double squared(double4 a) { return dot(a, a); }

inline double2 vlength(double2 a) { return sqrt(vdot(a, a)); }
inline double4 vlength(double4 a) { return sqrt(vdot(a, a)); }
inline double length(double2 a) { return vlength(a).x; }
inline double length(double4 a) { return vlength(a).x; }

// TODO div(sqrt) might be slower than mul(rsqrt)
inline auto normalize(double2 a) { return a / vlength(a); }
inline auto normalize(double4 a) { return a / vlength(a); }

template<typename V> bool is_unit(V v) { return ordered(1-1e-12, squared(v), 1+1e-12); }

template<int M> double2 round(double2 v, int mode) { return _mm_round_pd(v, M | _MM_FROUND_NO_EXC); }
template<int M> double4 round(double4 v, int mode) { return _mm256_round_pd(v, M | _MM_FROUND_NO_EXC); }

template<typename V> V round_nearest(V v) { return round<_MM_FROUND_TO_NEAREST_INT>(v); }
template<typename V> V round_zero(V v) { return round<_MM_FROUND_TO_ZERO>(v); }
template<typename V> V round_up(V v) { return round<_MM_FROUND_TO_POS_INF>(v); }
template<typename V> V round_down(V v) { return round<_MM_FROUND_TO_NEG_INF>(v); }

constexpr ulong SignBit64 = ulong(1) << 63;
inline double2 abs(double2 v) { return bit_and(v, lbroad2(~SignBit64)); }
inline double4 abs(double4 v) { return bit_and(v, lbroad4(~SignBit64)); }

inline double4 any_normal(double4 v) {
	double4 a = abs(v);
	if (a.x <= a.y && a.x <= a.z)
		return {0, -v.z, v.y, 0};
	if (a.y <= a.z)
		return {-v.z, 0, v.x, 0};
	return {-v.y, v.x, 0, 0};
}

inline double cross(double2 a, double2 b) { return a.x * b.y - b.x * a.y; }

inline double4 cross(double4 a, double4 b) {
	assert(a.w == 0);
	assert(b.w == 0);
	return {a.y * b.z - b.y * a.z,
            a.z * b.x - b.z * a.x,
            a.x * b.y - b.x * a.y,
			0};
}

// returns angle in range [0, PI)
inline auto angle(double4 a, double4 b) {
	return std::atan2(length(cross(a, b)), dot(a, b));
}

// solid angle between triangle and origin
inline double solid_angle(double4 A, double4 B, double4 C) {
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
inline double det(double4 a, double4 b, double4 c) {
	return a.x * det(b.yz, c.yz)
         - b.x * det(a.yz, c.yz)
		 + c.x * det(a.yz, b.yz);
}

// det(A) == det(transposed(A))
inline double det(double4 a, double4 b, double4 c, double4 d) {
	return a.x * det(b.yzww, c.yzww, d.yzww)
         - b.x * det(a.xzww, c.xzww, d.xzww)
		 + c.x * det(a.xyww, b.xyww, d.xyww)
		 - d.x * det(a.xyzw, b.xyzw, c.xyzw);
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
double4 uniform3(RND& rnd, double min, double max) {
	return {uniform(rnd, min, max), uniform(rnd, min, max), uniform(rnd, min, max), 1};
}

template<typename RND>
double4 uniform3v(RND& rnd, double min, double max) {
	return {uniform(rnd, min, max), uniform(rnd, min, max), uniform(rnd, min, max), 0};
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
double4 normal3p(RND& rnd, double mean, double stdev) {
	return {normal(rnd, mean, stdev), normal(rnd, mean, stdev), normal(rnd, mean, stdev), 1};
}
template<typename RND>
double4 normal3v(RND& rnd, double mean, double stdev) {
	return {normal(rnd, mean, stdev), normal(rnd, mean, stdev), normal(rnd, mean, stdev), 0};
}

template<typename RND>
double4 normal4(RND& rnd, double mean, double stdev) {
	return {normal(rnd, mean, stdev), normal(rnd, mean, stdev), normal(rnd, mean, stdev), normal(rnd, mean, stdev)};
}

// uniform on the unit sphere surface
template<typename RND>
double4 uniform_dir3(RND& rnd) { return normalize(normal3v(rnd, 0, 1)); }
template<typename RND>
double4 uniform_dir4(RND& rnd) { return normalize(normal4(rnd, 0, 1)); }

inline bool colinear(double4 a, double4 b, double4 c) {
	return squared(cross(b - a, c - a)) <= squared(1e-6);
}

// TODO define struct vec3
// w can always be assumed to be 0
// with alignas(32)
// operator ==
// hash

/*struct alignas(32) vec3 {
private:
	double4 s;
};

struct alignas(32) point3 {
	point3() : s{NAN, NAN, NAN, 1} { }
	explicit point3(double4 v) : s{v.x, v.y, v.z, 1} { }
	static point3 from_vec(double4 v) { return v; }
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
	point3 operator+(double4 v) const { return s + double4{v.x, v.y, v.z, 0}; }
	point3 operator-(double4 v) const { return s - double4{v.x, v.y, v.z, 0}; }

	point3& operator+=(double4 v) { s += v; return *this; }
	point3& operator-=(double4 v) { s -= v; return *this; }
	point3& operator+=(double4 v) { s += {v.x, v.y, v.z, 0}; return *this; }
	point3& operator-=(double4 v) { s -= {v.x, v.y, v.z, 0}; return *this; }

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

	double4 vec() const { return s; }
private:
	point3(double4 v) : s(v) { }
	double4 s;
};*/

/*inline void format_e(string& s, string_view spec, point3 v) {
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

}*/

using point3 = double4;

inline point3 avg(point3 a, point3 b) { return (a + b) / 2; }
inline float dot(point3 a, float4 b) { return dot(vconvert(a, float4), b); }
inline float dot(float4 a, point3 b) { return dot(a, vconvert(b, float4)); }

inline double4 compute_normal(point3 a, point3 b, point3 c) { return cross(b - a, c - a); }

inline void remove_dups(aligned_vector<double4>& vertices) {
	::remove_dups(vertices,
		static_cast<bool(*)(double4, double4)>(&lex_less),
		static_cast<bool(*)(double4, double4)>(&equal));
}

inline double4 point(double x, double y, double z) { return {x, y, z, 1}; }
