#pragma once
#include "common.h"
#include "format.h"
#include "hash.h"
#include <immintrin.h>

#define vshuffle __builtin_shufflevector
#define vconvert __builtin_convertvector

// use double4 for computation (with w=1 for points and w=0 for directions)

// 8x vec3 in component order,
// each component fits into one 256-bit AVX2 register
struct vec3_8 {
    float8 x, y, z;
};

struct fmat34 {
    float4 a, b, c;
};

// TODO rename to transform
struct mat34 {
    double4 a, b, c;
};

template<typename T>
struct vec_info {
};

#define VEC_INFO(T, N) template<> struct vec_info<T##N> { static int dim() { return N; } using Type = T; };

#define EQUAL2(T) inline bool equal(T a, T b) { return a.x == b.x && a.y == b.y; }
#define EQUAL3(T) inline bool equal(T a, T b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
#define EQUAL4(T) inline bool equal(T a, T b) { return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w; }

// TODO consider using crc32 avx instruction for hash code
#define STD_HASH(X) namespace std { template<> struct hash X; }
#define HASH2(T) STD_HASH(<T##2> { size_t operator()(T##2 v) const { return ::hash(v.x, v.y); } })
#define HASH3(T) STD_HASH(<T##3> { size_t operator()(T##3 v) const { return ::hash(v.x, v.y, v.z); } })
#define HASH4(T) STD_HASH(<T##4> { size_t operator()(T##4 v) const { return ::hash(v.x, v.y, v.z, v.w); } })

#define VEC(T, N) EQUAL##N(T##N); VEC_INFO(T, N); HASH##N(T);
#define VECX(N) VEC(int, N); VEC(long, N); VEC(cent, N); VEC(float, N); VEC(double, N);

VECX(2);
VECX(3);
VECX(4);

template<typename T>
struct equal_t {
	bool operator()(T a, T b) const {
		return equal(a, b);
	}
};

// SIMD helpers
inline double3 d3(double4 v) { return v.xyz; }
inline double4 d4(double3 v) { return vshuffle(v, v, 0, 1, 2, -1); }

inline int8 vmin(int8 a, int8 b) { return _mm256_min_epi32(a, b); }
inline int4 vmin(int4 a, int4 b) { return _mm_min_epi32(a, b); }

inline float8 vmin(float8 a, float8 b) { return _mm256_min_ps(a, b); }
inline float4 vmin(float4 a, float4 b) { return _mm_min_ps(a, b); }

inline double2 vmin(double2 a, double2 b) { return _mm_min_pd(a, b); }
inline double3 vmin(double3 a, double3 b) { return d3(_mm256_min_pd(d4(a), d4(b))); }
inline double4 vmin(double4 a, double4 b) { return _mm256_min_pd(a, b); }

inline int8 vmax(int8 a, int8 b) { return _mm256_max_epi32(a, b); }
inline int4 vmax(int4 a, int4 b) { return _mm_max_epi32(a, b); }

inline float8 vmax(float8 a, float8 b) { return _mm256_max_ps(a, b); }
inline float4 vmax(float4 a, float4 b) { return _mm_max_ps(a, b); }

inline double2 vmax(double2 a, double2 b) { return _mm_max_pd(a, b); }
inline double3 vmax(double3 a, double3 b) { return d3(_mm256_max_pd(d4(a), d4(b))); }
inline double4 vmax(double4 a, double4 b) { return _mm256_max_pd(a, b); }

inline int hmin(int8 a) {
     auto b = vmin(a, vshuffle(a, a, 2, 3, -1, -1, 6, 7, -1, -1));
     auto c = vmin(b, vshuffle(b, b, 1, -1, -1, -1, 5, -1, -1, -1));
     auto d = vmin(c, vshuffle(c, c, 4, -1, -1, -1, -1, -1, -1, -1));
     return d.x;
}

inline float hmin(float8 a) {
     auto b = vmin(a, vshuffle(a, a, 2, 3, -1, -1, 6, 7, -1, -1));
     auto c = vmin(b, vshuffle(b, b, 1, -1, -1, -1, 5, -1, -1, -1));
     auto d = vmin(c, vshuffle(c, c, 4, -1, -1, -1, -1, -1, -1, -1));
     return d.x;
}

inline int hmax(int8 a) {
     auto b = vmax(a, vshuffle(a, a, 2, 3, -1, -1, 6, 7, -1, -1));
     auto c = vmax(b, vshuffle(b, b, 1, -1, -1, -1, 5, -1, -1, -1));
     auto d = vmax(c, vshuffle(c, c, 4, -1, -1, -1, -1, -1, -1, -1));
     return d.x;
}

inline float hmax(float8 a) {
     auto b = vmax(a, vshuffle(a, a, 2, 3, -1, -1, 6, 7, -1, -1));
     auto c = vmax(b, vshuffle(b, b, 1, -1, -1, -1, 5, -1, -1, -1));
     auto d = vmax(c, vshuffle(c, c, 4, -1, -1, -1, -1, -1, -1, -1));
     return d.x;
}

inline int4 as_int(float4 a) { return _mm_cvtps_epi32(a); }
inline float4 as_float(int4 a) { return _mm_cvtepi32_ps(a); }

// Note: destination must be one of the arguments, otherwise it will not compile to 1 instruction
inline float4 fma(float4 a, float4 b, float4 c) { return _mm_fmadd_ps(a, b, c); }
inline float8 fma(float8 a, float8 b, float8 c) { return _mm256_fmadd_ps(a, b, c); }
inline double4 fma(double4 a, double4 b, double4 c) { return _mm256_fmadd_pd(a, b, c); }

inline bool all(int2 v) { return v[0] && v[1]; }
inline bool all(int3 v) { return v[0] && v[1] && v[2]; }
inline bool all(int4 v) { return v[0] && v[1] && v[2] && v[3]; }

inline bool any(int2 v) { return v[0] || v[1]; }
inline bool any(int3 v) { return v[0] || v[1] || v[2]; }
inline bool any(int4 v) { return v[0] || v[1] || v[2] || v[3]; }
