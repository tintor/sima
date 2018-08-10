#pragma once
#include "glm.h"
#include "scalar.h"
#include <random>

template<typename RandomEngine>
ivec3 random_vector(RandomEngine& rnd, int a, int b) {
	std::uniform_int_distribution<int> dist(a, b);
	return ivec3{dist(rnd), dist(rnd), dist(rnd)};
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
inline llvec3 muli(llvec3 a, cent b) {
	return {muli(a.x, b), muli(a.y, b), muli(a.z, b)};
}

// throws std::overflow_error in case of int64 overflow
inline long doti(ivec3 a, ivec3 b) {
	long x = (long)a.x * (long)b.x;
	long y = (long)a.y * (long)b.y;
	long z = (long)a.z * (long)b.z;
	return addi(addi(x, y), z);
}
inline long doti(lvec3 a, lvec3 b) {
	long x = muli(a.x, b.x);
	long y = muli(a.y, b.y);
	long z = muli(a.z, b.z);
	return addi(addi(x, y), z);
}
inline long doti(ivec3 a, lvec3 b) {
	return doti(vconvert(a, lvec3), b);
}
inline long doti(lvec3 a, ivec3 b) {
	return doti(a, vconvert(b, lvec3));
}
inline cent doti(llvec3 a, llvec3 b) {
	cent x = muli(a.x, b.x);
	cent y = muli(a.y, b.y);
	cent z = muli(a.z, b.z);
	return addi(addi(x, y), z);
}

// throws std::overflow_error in case of int64 overflow
inline long squaredi(lvec3 a) {
	return doti(a, a);
}
inline cent squaredi(llvec3 a) {
	return doti(a, a);
}

// throws std::overflow_error in case of int64 overflow
inline lvec3 crossi(lvec3 a, lvec3 b) {
	long x = subi(muli(a.y, b.z), muli(a.z, b.y));
	long y = subi(muli(a.z, b.x), muli(a.x, b.z));
	long z = subi(muli(a.x, b.y), muli(a.y, b.x));
	return {x, y, z};
}
inline llvec3 crossi(llvec3 a, llvec3 b) {
	cent x = subi(muli(a.y, b.z), muli(a.z, b.y));
	cent y = subi(muli(a.z, b.x), muli(a.x, b.z));
	cent z = subi(muli(a.x, b.y), muli(a.y, b.x));
	return {x, y, z};
}
inline lvec3 crossi(ivec3 a, ivec3 b) {
	return crossi(vconvert(a, lvec3), vconvert(b, lvec3));
}

// throws std::overflow_error in case of int64 overflow
inline lvec3 normali(ivec3 a, ivec3 b, ivec3 c) {
	return crossi(subi(b, a), subi(c, a));
}
inline llvec3 normali(lvec3 a, lvec3 b, lvec3 c) {
	llvec3 A = vconvert(a, llvec3);
	llvec3 B = vconvert(b, llvec3);
	llvec3 C = vconvert(c, llvec3);
	return crossi(subi(B, A), subi(C, A));
}

// throws std::overflow_error in case of int64 overflow
inline bool colinear(ivec3 a, ivec3 b, ivec3 c) {
	return equal(crossi(subi(b, a), subi(c, a)), lvec3{0, 0, 0});
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

inline lvec3 any_normal(lvec3 v) {
	// TODO if z == long::min and y != long::min then negi(y) instead
	if (v.x == 0)
		return {0, negi(v.z), v.y};
	if (v.y == 0)
		return {negi(v.z), 0, v.x};
	if (v.z == 0)
		return {negi(v.y), v.x, 0};
	THROW(invalid_argument);
}
