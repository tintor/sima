#pragma once
#include "common.h"
#include "format.h"

#define vshuffle __builtin_shufflevector
#define vconvert __builtin_convertvector

using ivec2 = int2;

using ivec3 = int3;
using lvec3 = long3;
using llvec3 = cent3;

using vec2 = float2;
using dvec2 = double2;

using vec3 = float3;
using dvec3 = double3;

using vec4 = float4;
using dvec4 = double4;

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

inline bool equal(ivec2 a, ivec2 b) {
	return a.x == b.x && a.y == b.y;
}

inline bool equal(ivec3 a, ivec3 b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool equal(lvec3 a, lvec3 b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool equal(llvec3 a, llvec3 b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool equal(vec3 a, vec3 b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool equal(dvec3 a, dvec3 b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

template<typename T>
struct equal_t {
	bool operator()(T a, T b) const {
		return equal(a, b);
	}
};

namespace std {

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

}

