#pragma once
#include "common.h"
#include "format.h"
#include <immintrin.h>

#define vshuffle __builtin_shufflevector
#define vconvert __builtin_convertvector

// use int4 for meshes (with w=1)
// use double4 for computation (with w=1 for points and w=0 for directions)

using ivec2 = int2;

using ivec3 = int3;
using lvec3 = long3;
using llvec3 = cent3;

using ivec4 = int4;

using vec2 = float2;
using dvec2 = double2;

using vec3 = float3;
using dvec3 = double3;

using vec4 = float4;
using dvec4 = double4;

using ivec8 = int8;
using vec8 = float8;

// 8x vec3 in component order,
// each component fits into one 256-bit AVX2 register
struct ivec3_8 {
    int8 x, y, z;
};

struct vec3_8 {
    float8 x, y, z;
};

struct mat34 {
    float4 a, b, c;
};

struct dmat34 {
    double4 a, b, c;
};

template<typename T>
struct vec_info {
};

template<>
struct vec_info<vec2> {
	static int dim() { return 2; }
	using Type = float;
};

template<>
struct vec_info<dvec2> {
	static int dim() { return 2; }
	using Type = double;
};

template<>
struct vec_info<ivec2> {
	static int dim() { return 2; }
	using Type = int;
};

template<>
struct vec_info<vec3> {
	static int dim() { return 3; }
	using Type = float;
};

template<>
struct vec_info<dvec3> {
	static int dim() { return 3; }
	using Type = double;
};

template<>
struct vec_info<ivec3> {
	static int dim() { return 3; }
	using Type = int;
};

#define EQUAL2(T) inline bool equal(T a, T b) { return a.x == b.x && a.y == b.y; }
#define EQUAL3(T) inline bool equal(T a, T b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
#define EQUAL4(T) inline bool equal(T a, T b) { return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w; }

EQUAL2(ivec2);

EQUAL3(ivec3);
EQUAL3(lvec3);
EQUAL3(llvec3);
EQUAL3(vec3);
EQUAL3(dvec3);

EQUAL4(ivec4);

template<typename T>
struct equal_t {
	bool operator()(T a, T b) const {
		return equal(a, b);
	}
};

namespace std {

// TODO consider using crc32 avx instruction for hash code
template<>
struct hash<ivec2> {
	size_t operator()(ivec2 v) const {
		size_t seed = 0;
		hash_combine(seed, v.x);
		hash_combine(seed, v.y);
		return seed;
	}
};

template<>
struct hash<ivec3> {
	size_t operator()(ivec3 v) const {
		size_t seed = 0;
		hash_combine(seed, v.x);
		hash_combine(seed, v.y);
		hash_combine(seed, v.z);
		return seed;
	}
};

template<>
struct hash<dvec3> {
	size_t operator()(dvec3 v) const {
		size_t seed = 0;
		hash_combine(seed, v.x);
		hash_combine(seed, v.y);
		hash_combine(seed, v.z);
		return seed;
	}
};

}

inline int8 vmin(int8 a, int8 b) { return _mm256_min_epi32(a, b); }
inline int4 vmin(int4 a, int4 b) { return _mm_min_epi32(a, b); }
inline float8 vmin(float8 a, float8 b) { return _mm256_min_ps(a, b); }
inline float4 vmin(float4 a, float4 b) { return _mm_min_ps(a, b); }

inline int8 vmax(int8 a, int8 b) { return _mm256_max_epi32(a, b); }
inline int4 vmax(int4 a, int4 b) { return _mm_max_epi32(a, b); }
inline float8 vmax(float8 a, float8 b) { return _mm256_max_ps(a, b); }
inline float4 vmax(float4 a, float4 b) { return _mm_max_ps(a, b); }

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
