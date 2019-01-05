#pragma once
#include <core/int.h>
#include <core/format.h>
#include <core/hash.h>
#include <immintrin.h>

#define vshuffle __builtin_shufflevector
#define vconvert __builtin_convertvector

#define m128 __m128d
#define m256 __m256d
#define m512 __m512d

inline int4 cast4(int3 a) { return vshuffle(a, a, 0, 1, 2, -1); }
inline long4 cast4(long3 a) { return vshuffle(a, a, 0, 1, 2, -1); }
inline float4 cast4(float3 a) { return vshuffle(a, a, 0, 1, 2, -1); }
inline double4 cast4(double3 a) { return vshuffle(a, a, 0, 1, 2, -1); }

// use m128 and m256 to allow input of differrent types (but equal size)
inline m128 bit_and(m128 a, m128 b) { return _mm_and_pd(a, b); }
inline m256 bit_and(m256 a, m256 b) { return _mm256_and_pd(a, b); }

inline m128 bit_or(m128 a, m128 b) { return _mm_or_pd(a, b); }
inline m256 bit_or(m256 a, m256 b) { return _mm256_or_pd(a, b); }

inline m128 bit_xor(m128 a, m128 b) { return _mm_xor_pd(a, b); }
inline m256 bit_xor(m256 a, m256 b) { return _mm256_xor_pd(a, b); }

inline bool all(int2 v) { return v.x != 0 && v.y != 0; }
inline bool all(int3 v) { return _mm_movemask_ps(cast4(v)) == 0xF; }
inline bool all(int4 v) { return _mm_movemask_ps(v) == 0xF; }

inline bool all(long2 v) { return _mm_movemask_pd(v) == 0x3; }
inline bool all(long3 v) { return _mm256_movemask_pd(cast4(v)) == 0xF; }
inline bool all(long4 v) { return _mm256_movemask_pd(v) == 0xF; }

inline bool any(int2 v) { return v.x != 0 || v.y != 0; }
inline bool any(int3 v) { return _mm_movemask_ps(cast4(v)) != 0; }
inline bool any(int4 v) { return _mm_movemask_ps(v) != 0; }

inline bool any(long2 v) { return _mm_movemask_pd(v) != 0; }
inline bool any(long3 v) { return _mm256_movemask_pd(cast4(v)) != 0; }
inline bool any(long4 v) { return _mm256_movemask_pd(v) != 0; }

inline bool equal(int2 a, int2 b) { return all(a == b); }
inline bool equal(int3 a, int3 b) { return all(a == b); }
inline bool equal(int4 a, int4 b) { return all(a == b); }

// requires AVX512
//inline bool equal(double2 a, double2 b) { return _mm_cmpeq_epu64_mask(a, b) == 0x3; }
//inline bool equal(double4 a, double4 b) { return _mm256_cmpeq_epu64_mask(a, b) == 0xF; }

inline bool equal(double2 a, double2 b) { return all(a == b); }
inline bool equal(double3 a, double3 b) { return all(a == b); }
inline bool equal(double4 a, double4 b) { return all(a == b); }

inline double2 floor(double2 a) { return _mm_floor_pd(a); }
inline double3 floor(double3 a) { double4 v = _mm256_floor_pd(cast4(a)); return v.xyz; }
inline double4 floor(double4 a) { return _mm256_floor_pd(a); }

inline double2 ceil(double2 a) { return _mm_ceil_pd(a); }
inline double3 ceil(double3 a) { double4 v = _mm256_ceil_pd(cast4(a)); return v.xyz; }
inline double4 ceil(double4 a) { return _mm256_ceil_pd(a); }

inline hash operator<<(hash h, int2 a) { return h << a.x << a.y; }

inline hash operator<<(hash h, double2 a) { return h << a.x << a.y; }
inline hash operator<<(hash h, double3 a) { return h << a.x << a.y << a.z; }
inline hash operator<<(hash h, double4 a) { return h << a.x << a.y << a.z << a.w; }

template<typename T>
struct equal_t {
	bool operator()(T a, T b) const {
		return equal(a, b);
	}
};

inline int8 vmin(int8 a, int8 b) { return _mm256_min_epi32(a, b); }
inline int4 vmin(int4 a, int4 b) { return _mm_min_epi32(a, b); }

inline float8 vmin(float8 a, float8 b) { return _mm256_min_ps(a, b); }
inline float4 vmin(float4 a, float4 b) { return _mm_min_ps(a, b); }

inline double2 vmin(double2 a, double2 b) { return _mm_min_pd(a, b); }
inline double3 vmin(double3 a, double3 b) { double4 e = _mm256_min_pd(cast4(a), cast4(b)); return e.xyz; }
inline double4 vmin(double4 a, double4 b) { return _mm256_min_pd(a, b); }

inline int8 vmax(int8 a, int8 b) { return _mm256_max_epi32(a, b); }
inline int4 vmax(int4 a, int4 b) { return _mm_max_epi32(a, b); }

inline float8 vmax(float8 a, float8 b) { return _mm256_max_ps(a, b); }
inline float4 vmax(float4 a, float4 b) { return _mm_max_ps(a, b); }

inline double2 vmax(double2 a, double2 b) { return _mm_max_pd(a, b); }
inline double3 vmax(double3 a, double3 b) { double4 e = _mm256_max_pd(cast4(a), cast4(b)); return e.xyz; }
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
