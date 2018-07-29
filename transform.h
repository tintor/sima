#pragma once

#include "glm.h"
#include <immintrin.h>

typedef float float3 __attribute__((ext_vector_type(3)));
typedef float float4 __attribute__((ext_vector_type(4)));
typedef int int4 __attribute__((ext_vector_type(4)));

// __builtin_convertvector

// converts float4 to int4
// __m128i _mm_cvtps_epi32 (__m128 a)
// convert int4 to float4
// __m128 _mm_cvtepi32_ps (__m128i a)

// transform vector by matrix with 5 instructions!
// +2 instructions if input/output is int4
inline float4 transform(float4 v, float4 m[3]) {
	// dot products
	float4 x = _mm_dp_ps(v, m[0], 0xF1);
	float4 y = _mm_dp_ps(v, m[1], 0xF2);
	float4 z = _mm_dp_ps(v, m[2], 0xF4);
	return x + y + z;
}

void translate(imesh3& m, float3 t);
void transform(imesh3& m, float4 t[3]);

struct transform3 {
	/*dmat3 orientation;
	dvec3 position;

    transform3(const dvec3& position, dquat orientation)
        : orientation(orientation), position(position) { }

    dvec3 to_global(dvec3 v) const { return orientation * v + position; }
    dvec3 to_local(dvec3 v) const { return (v - position) * orientation; }

    segment3 to_local(const segment3& p) const {
		return segment3(to_local(p.a), to_local(p.b));
	}

	dvec3 to_global_dir(const dvec3& p) const {
		return orientation * p;
	}

	segment3 to_global(const segment3& p) const {
		return segment3(to_global(p.a), to_global(p.b));
	}

	triangle3 to_local(const triangle3& p) const {
		return triangle3(to_local(p.a), to_local(p.b), to_local(p.c));
	}

	triangle3 to_global(const triangle3& p) const {
		return triangle3(to_global(p.a), to_global(p.b), to_global(p.c));
	}*/
};

transform3 combine(const transform3& a, const transform3& b);

/*
transform3 combine(const transform3& a, const transform3& b) {
    // TODO create two dmat4 from two transform3
	dmat4 a4, b4;
    // TODO figure out combination formula: ie: a * inverse(b)
	dmat4 c4 = a4 * inverse(b4);
    // TODO convert c4 back to transform3
	throw new std::runtime_error("combine() unimplemented");
}
*/
