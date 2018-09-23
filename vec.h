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

inline bool all(int2 v) { return v.x != 0 && v.y != 0; }
inline bool all(int4 v) { return _mm_movemask_ps(v) == 0xF; }
inline bool any(int2 v) { return v.x != 0 || v.y != 0; }
inline bool any(int4 v) { return _mm_movemask_ps(v) != 0; }

inline bool all(long2 v) { return _mm_movemask_pd(v) == 0x3; }
inline bool all(long4 v) { return _mm256_movemask_pd(v) == 0xF; }
inline bool any(long2 v) { return _mm_movemask_pd(v) != 0; }
inline bool any(long4 v) { return _mm256_movemask_pd(v) != 0; }

inline bool equal(int2 a, int2 b) { return all(a == b); }
inline bool equal(double2 a, double2 b) { return all(a == b); }
inline bool equal(double4 a, double4 b) { return all(a == b); }

inline double4 floor(double4 a) { return _mm256_floor_pd(a); }
inline double4 ceil(double4 a) { return _mm256_ceil_pd(a); }

// requires AVX512
//inline bool equal(double2 a, double2 b) { return _mm_cmpeq_epu64_mask(a, b) == 0x3; }
//inline bool equal(double4 a, double4 b) { return _mm256_cmpeq_epu64_mask(a, b) == 0xF; }

inline hash operator<<(hash h, int2 a) { return h << a.x << a.y; }
inline hash operator<<(hash h, double2 a) { return h << a.x << a.y; }
inline hash operator<<(hash h, double4 a) { return h << a.x << a.y << a.z << a.w; }

template<typename T>
struct equal_t {
	bool operator()(T a, T b) const {
		return equal(a, b);
	}
};

// SIMD helpers
inline int8 vmin(int8 a, int8 b) { return _mm256_min_epi32(a, b); }
inline int4 vmin(int4 a, int4 b) { return _mm_min_epi32(a, b); }

inline float8 vmin(float8 a, float8 b) { return _mm256_min_ps(a, b); }
inline float4 vmin(float4 a, float4 b) { return _mm_min_ps(a, b); }

inline double2 vmin(double2 a, double2 b) { return _mm_min_pd(a, b); }
inline double4 vmin(double4 a, double4 b) { return _mm256_min_pd(a, b); }

inline int8 vmax(int8 a, int8 b) { return _mm256_max_epi32(a, b); }
inline int4 vmax(int4 a, int4 b) { return _mm_max_epi32(a, b); }

inline float8 vmax(float8 a, float8 b) { return _mm256_max_ps(a, b); }
inline float4 vmax(float4 a, float4 b) { return _mm_max_ps(a, b); }

inline double2 vmax(double2 a, double2 b) { return _mm_max_pd(a, b); }
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
