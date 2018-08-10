#pragma once
#include "int.h"
#include "exception.h"
#include "bits_util.h"

inline int sign(int a) {
	return (0 < a) - (a < 0);
}

inline int sign(long a) {
	return (0 < a) - (a < 0);
}

inline int sign(int128 a) {
	return (0 < a) - (a < 0);
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
