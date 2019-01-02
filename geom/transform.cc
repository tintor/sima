#include <geom/transform.h>
#include <core/format.h>
#include <core/util.h>
#include <core/exception.h>

// 7 instructions
// input a.w is ignored (as long as matrix is fully initialized 4x4)
// output w is 0
float4 transform(float4 a, float4 m[4]) {
    float4 b = m[3];
    b = fma(a.xxxx, m[0], b);
    b = fma(a.yyyy, m[1], b);
    b = fma(a.zzzz, m[2], b);
    return b;
}

// transform 2 vectors at once
// 7 instructions!
// 3.5 instructions per vector
// m must be pre-expanded
float8 transform(float8 a, float8 m[4]) {
    float8 b = m[3];
    b = fma(vshuffle(a, a, 0, 0, 0, 0, 4, 4, 4, 4), m[0], b);
    b = fma(vshuffle(a, a, 1, 1, 1, 1, 5, 5, 5, 5), m[1], b);
    b = fma(vshuffle(a, a, 2, 2, 2, 2, 6, 6, 6, 6), m[2], b);
    return b;
}

// inner loop is 15 instructions for 8 vectors!
/*void transform(cspan<vec3_8> in, fmat34 m, span<vec3_8> out) {
	if (in.size() != out.size())
		THROW(invalid_argument, "array sizes must be the same");
	if (!aligned<32>(in.begin()) || !aligned<32>(out.begin()))
		THROW(invalid_argument, "arrays must be 32-byte aligned");

	float8 max = m.a.xxxxxxxx;
    float8 may = m.a.yyyyyyyy;
    float8 maz = m.a.zzzzzzzz;
    float8 maw = m.a.wwwwwwww;

    float8 mbx = m.b.xxxxxxxx;
    float8 mby = m.b.yyyyyyyy;
    float8 mbz = m.b.zzzzzzzz;
    float8 mbw = m.b.wwwwwwww;

    float8 mcx = m.c.xxxxxxxx;
    float8 mcy = m.c.yyyyyyyy;
    float8 mcz = m.c.zzzzzzzz;
    float8 mcw = m.c.wwwwwwww;

	auto a = in.begin();
	auto b = out.begin();
	while (a < in.end()) {
	    float8 v = maw;
	    v = fma(a->x, max, v);
	    v = fma(a->y, may, v);
	    v = fma(a->z, maz, v);
	    b->x = v;

	    v = mbw;
   		v = fma(a->x, mbx, v);
	    v = fma(a->y, mby, v);
	    v = fma(a->z, mbz, v);
	    b->y = v;

	    v = mcw;
	    v = fma(a->x, mcx, v);
	    v = fma(a->y, mcy, v);
	    v = fma(a->z, mcz, v);
	    b->z = v;

		a += 1;
		b += 1;
	}
}*/

/*void translate(span<vec3_8> vectors, ivec4 t) {
	if (!aligned<32>(vectors.begin()))
		THROW(invalid_argument, "address must be 32-byte aligned");

	ivec8 tx = t.xxxxxxxx;
	ivec8 ty = t.yyyyyyyy;
	ivec8 tz = t.zzzzzzzz;

	for(auto it = vectors.begin(); it < vectors.end(); it++) {
		it->x += tx;
		it->y += ty;
		it->z += tz;
	}
}

void translate(span<ivec4> vectors, ivec4 offset) {
	ivec8* it = reinterpret_cast<ivec8*>(vectors.begin());
	ivec8* end = it + vectors.size() / 2;
	if (!aligned<32>(vectors.begin()))
		THROW(invalid_argument, "address must be 32-byte aligned");

	ivec8 offset2 = offset.xyzwxyzw;
	while (it < end)
		*it++ += offset2;

	if (vectors.size() % 2 == 1)
		vectors.back() += offset;
}*/
