#pragma once
#include "triangle.h"
#include "array_ptr.h"

inline float4 transform(float4 v, mat34 matrix) {
	// dot products
	float4 x = _mm_dp_ps(v, matrix.a, 0b11110001);
	float4 y = _mm_dp_ps(v, matrix.b, 0b11110010);
	float4 z = _mm_dp_ps(v, matrix.c, 0b11110100);
	float4 q = x + y + z;
	q.w = v.w;
	return q;
}

void transform(array_cptr<vec3_8> in, mat34 matrix, array_ptr<vec3_8> out);

void translate(array_ptr<ivec4> vectors, ivec4 t);
void translate(array_ptr<ivec3_8> vectors, ivec4 t);

/*inline void translate(imesh3& mesh, ivec3 t) {
	ivec3* begin = &mesh.begin()->a;
	ivec3* end = &mesh.begin()->a;
	translate(array_ptr<ivec3>(begin, end), t);
}*/

//void transform(imesh3& mesh, float4 t[3]);

struct transform3 {
	mat34 matrix;
	/*

    transform3(const dvec4& position, dvec4 orientation)
        : orientation(orientation), position(position) { }

    dvec3 to_global(dvec3 v) const { return orientation * v + position; }
    dvec3 to_local(dvec3 v) const { return (v - position) * orientation; }

    isegment3 to_local(const isegment3& p) const {
		return isegment3(to_local(p.a), to_local(p.b));
	}

	dvec3 to_global_dir(const dvec3& p) const {
		return orientation * p;
	}

	isegment3 to_global(const isegment3& p) const {
		return isegment3(to_global(p.a), to_global(p.b));
	}

	itriangle3 to_local(const itriangle3& p) const {
		return itriangle3(to_local(p.a), to_local(p.b), to_local(p.c));
	}

	itriangle3 to_global(const itriangle3& p) const {
		return itriangle3(to_global(p.a), to_global(p.b), to_global(p.c));
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
