#pragma once

#define GLM_FORCE_RADIANS
#define GLM_FORCE_AVX2
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

#include "common.h"

using glm::ivec2;
using glm::ivec3;
using glm::dvec3;
using lvec3 = glm::tvec3<long>;
using llvec3 = glm::tvec3<int128>;
using glm::dmat3;

#include "range.h"
#include "bits.h"
#include "exception.h"
#include <ostream>

namespace glm {

inline void format_e(std::string& s, std::string_view spec, ivec2 v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
}

inline void format_e(std::string& s, std::string_view spec, ivec3 v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
	s += ' ';
	::format_e(s, "", v.z);
}

inline void format_e(std::string& s, std::string_view spec, lvec3 v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
	s += ' ';
	::format_e(s, "", v.z);
}

inline void format_e(std::string& s, std::string_view spec, llvec3 v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
	s += ' ';
	::format_e(s, "", v.z);
}

}

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

template<typename T>
T gcd_unsigned(T u, T v) {
    if (u == 0) return v;
    if (v == 0) return u;
    int shift = ctz(u | v);
    u >>= ctz(u);
    do {
        v >>= ctz(v);
        if (u > v) {
            auto t = v;
            v = u;
            u = t;
        }
        v -= u;
    } while (v != 0);
    return u << shift;
}

inline long gcd(lvec3 s) {
	ulong x = std::abs(s.x);
	ulong y = std::abs(s.y);
	ulong z = std::abs(s.z);

	// TODO verify if sorting helps
	ulong a, b, c;
	if (x <= y) {
		if (z <= x) {
			a = z;
			b = x;
			c = y;
		} else if (y <= z) {
			a = x;
			b = y;
			c = z;
		} else {
			a = x;
			b = z;
			c = y;
		}
	} else {
		if (z <= y) {
			a = z;
			b = y;
			c = x;
		} else if (x <= z) {
			a = y;
			b = x;
			c = z;
		} else {
			a = y;
			b = z;
			c = x;
		}
	}
	// compute gdc of smallest numbers first
	return gcd_unsigned(c, gcd_unsigned(b, a));
}

inline long addi(int a, int b) {
	return (long)a + (long)b;
}
#define ADDI(R, A, B) inline R addi(A a, B b) { \
	R e; \
	if (__builtin_add_overflow(a, b, &e)) \
		THROW(overflow_error, "addi(%s, %s)", a, b); \
	return e; \
}
ADDI(long, long, int);
ADDI(long, int, long);
ADDI(long, long, long);
ADDI(int128, long, int128);
ADDI(int128, int128, long);
ADDI(int128, int128, int128);

inline long negi(int a) {
	return -(long)a;
}
inline long negi(long a) {
	if (a == std::numeric_limits<long>::min())
		THROW(overflow_error, "negi(%s)", a);
	return -a;
}
inline int128 negi(int128 a) {
	if (a == std::numeric_limits<int128>::min())
		THROW(overflow_error, "negi(%s)", a);
	return -a;
}

inline long subi(int a, int b) {
	return (long)a - (long)b;
}
#define SUBI(R, A, B) inline R subi(A a, B b) { \
	R e; \
	if (__builtin_sub_overflow(a, b, &e)) \
		THROW(overflow_error, "subi(%s, %s)", a, b); \
	return e; \
}
SUBI(long, int, long);
SUBI(long, long, int);
SUBI(long, long, long);
SUBI(int128, long, int128);
SUBI(int128, int128, long);
SUBI(int128, int128, int128);

inline long muli(int a, int b) {
	return (long)a * (long)b;
}
#define MULI(R, A, B) inline R muli(A a, B b) { \
	R e; \
	if (__builtin_mul_overflow(a, b, &e)) \
		THROW(overflow_error, "muli(%s, %s)", a, b); \
	return e; \
}
MULI(long, int, long);
MULI(long, long, int);
MULI(long, long, long);
MULI(int128, long, int128);
MULI(int128, int128, long);
MULI(int128, int128, int128);

inline lvec3 addi(ivec3 a, ivec3 b) {
	return {addi(a.x, b.x), addi(a.y, b.y), addi(a.z, b.z)};
}
inline lvec3 addi(lvec3 a, ivec3 b) {
	return {addi(a.x, b.x), addi(a.y, b.y), addi(a.z, b.z)};
}
inline lvec3 addi(ivec3 a, lvec3 b) {
	return {addi(a.x, b.x), addi(a.y, b.y), addi(a.z, b.z)};
}
inline lvec3 addi(lvec3 a, lvec3 b) {
	return {addi(a.x, b.x), addi(a.y, b.y), addi(a.z, b.z)};
}
inline llvec3 addi(llvec3 a, llvec3 b) {
	return {addi(a.x, b.x), addi(a.y, b.y), addi(a.z, b.z)};
}

inline lvec3 subi(ivec3 a, ivec3 b) {
	return {subi(a.x, b.x), subi(a.y, b.y), subi(a.z, b.z)};
}
inline lvec3 subi(lvec3 a, lvec3 b) {
	return {subi(a.x, b.x), subi(a.y, b.y), subi(a.z, b.z)};
}
inline llvec3 subi(llvec3 a, llvec3 b) {
	return {subi(a.x, b.x), subi(a.y, b.y), subi(a.z, b.z)};
}

inline lvec3 muli(lvec3 a, long b) {
	return {muli(a.x, b), muli(a.y, b), muli(a.z, b)};
}
inline llvec3 muli(llvec3 a, int128 b) {
	return {muli(a.x, b), muli(a.y, b), muli(a.z, b)};
}

// throws std::overflow_error in case of int64 overflow
inline long doti(lvec3 a, lvec3 b) {
	long x = muli(a.x, b.x);
	long y = muli(a.y, b.y);
	long z = muli(a.z, b.z);
	return addi(addi(x, y), z);
}
inline int128 doti(llvec3 a, llvec3 b) {
	int128 x = muli(a.x, b.x);
	int128 y = muli(a.y, b.y);
	int128 z = muli(a.z, b.z);
	return addi(addi(x, y), z);
}

// throws std::overflow_error in case of int64 overflow
inline long squaredi(lvec3 a) {
	return doti(a, a);
}
inline int128 squaredi(llvec3 a) {
	return doti(a, a);
}

// throws std::overflow_error in case of int64 overflow
inline lvec3 crossi(lvec3 a, lvec3 b) {
	long x = subi(muli(a.y, b.z), muli(a.z, b.y));
	long y = subi(muli(a.z, b.x), muli(a.x, b.z));
	long z = subi(muli(a.x, b.y), muli(a.y, b.x));
	return {x, y, z};
}

// throws std::overflow_error in case of int64 overflow
inline lvec3 normali(ivec3 a, ivec3 b, ivec3 c) {
	return crossi(subi(b, a), subi(c, a));
}

// throws std::overflow_error in case of int64 overflow
inline bool colinear(ivec3 a, ivec3 b, ivec3 c) {
	return crossi(subi(b, a), subi(c, a)) == lvec3(0, 0, 0);
}

template<typename Vec>
bool add_overflow(Vec a, Vec b, Vec& c) {
	return __builtin_add_overflow(a.x, b.x, &c.x)
		|| __builtin_add_overflow(a.y, b.y, &c.y)
		|| __builtin_add_overflow(a.z, b.z, &c.z);
}

template<typename Vec>
bool sub_overflow(Vec a, Vec b, Vec& c) {
	return __builtin_sub_overflow(a.x, b.x, &c.x)
		|| __builtin_sub_overflow(a.y, b.y, &c.y)
		|| __builtin_sub_overflow(a.z, b.z, &c.z);
}
