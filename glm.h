#pragma once

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

#include "common.h"
#include "range.h"
#include <ostream>

namespace glm {

template<typename T>
std::ostream& operator<<(std::ostream& os, const tvec3<T>& v) { return os << v.x << ' ' << v.y << ' ' << v.z; }

template<typename T>
std::ostream& operator<<(std::ostream& os, const tvec2<T>& v) { return os << v.x << ' ' << v.y; }

inline void format_e(std::string& s, std::string_view spec, ivec3 v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
	s += ' ';
	::format_e(s, "", v.z);
}

inline void format_e(std::string& s, std::string_view spec, ivec2 v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
}

}

using lvec3 = glm::tvec3<long>;

namespace std {

template<typename T>
struct hash<glm::tvec2<T>> {
	size_t operator()(const glm::tvec2<T>& v) const {
		size_t seed = 0;
		hash_combine(seed, v.x);
		hash_combine(seed, v.y);
		return seed;
	}
};

template<typename T>
struct hash<glm::tvec3<T>> {
	size_t operator()(const glm::tvec3<T>& v) const {
		size_t seed = 0;
		hash_combine(seed, v.x);
		hash_combine(seed, v.y);
		hash_combine(seed, v.z);
		return seed;
	}
};

}

inline long mul(int a, int b) {
	// TODO maybe faster with assembly
	return static_cast<long>(a) * static_cast<long>(b);
}

inline int128 mul(long a, long b) {
	// TODO maybe faster with assembly
	return static_cast<int128>(a) * static_cast<int128>(b);
/*    int128 res;
    __asm__ (
        "movq  %1, %%rax;\n\t"          // rax = a
        "movl  %3, %%ecx;\n\t"          // ecx = b
        "imulq %2;\n\t"                 // rdx:rax = a * b
        "shrdq %%cl, %%rdx, %%rax;\n\t" // rax = int64_t (rdx:rax >> s)
        "movq  %%rax, %0;\n\t"          // res = rax
        : "=rm" (res)
        : "rm"(a), "rm"(b)
        : "%rax", "%rdx", "%ecx");
    return res;*/
}

// TODO use asserts to verify no overflow
inline int128 dot(lvec3 a, lvec3 b) {
	return mul(a.x, b.x) + mul(a.y, b.y) + mul(a.z, b.z);
}

inline int128 squared(lvec3 a) {
	return dot(a, a);
}

// TODO use asserts to verify no overflow
// TODO use clang __builtin_sadd_overflow extensions
// takes int32_t vectors and returns int64_t vector (no overflow)
inline lvec3 cross(ivec3 a, ivec3 b) {
	auto x = mul(a.y, b.z) - mul(a.z, b.y);
	auto y = mul(a.z, b.x) - mul(a.x, b.z);
	auto z = mul(a.x, b.y) - mul(a.y, b.x);
	return lvec3(x, y, z);
}

