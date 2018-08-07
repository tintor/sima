#include "properties.h"
#include "aabb.h"
#include "scalar.h"
#include "ivec.h"
#include "primitives.h"

double triangle_volume(ivec3 a, ivec3 b, ivec3 c) {
	double z = a.z;
	z += b.z;
	z += c.z;
	double s;
	s = ((double)a.y + b.y) * ((double)a.x - b.x);
	s += ((double)b.y + c.y) * ((double)b.x - c.x);
	s += ((double)c.y + a.y) * ((double)c.x - a.x);
	return s * z;
}

double signed_volume(ivec3 a, ivec3 b, ivec3 c, ivec3 d) {
	double v;
	v = triangle_volume(a, b, c);
	v += triangle_volume(d, b, a);
	v += triangle_volume(a, c, b);
	v += triangle_volume(d, a, c);
	return v / 6.0;
}

double signed_volume(const imesh3& mesh) {
	double v = 0;
	for (const itriangle3& f : mesh)
		v += triangle_volume(f.a, f.b, f.c);
	return v / 6.0;
}

double volume(const imesh3& mesh) {
	return std::abs(signed_volume(mesh));
}

// TODO test with volume of cube (randomly rotated)

// Center of mass of a valid polyhedron
dvec3 center_of_mass(const imesh3& mesh) {
	dvec3 P = {0, 0, 0};
	double V = 0;
	for (const itriangle3& f : mesh) {
		dvec3 a = vconvert(f.a, double3);
		dvec3 b = vconvert(f.b, double3);
	   	dvec3 c = vconvert(f.c, double3);
		double v = dot(a, cross(b, c));
		P += (a + b + c) * v;
		V += v;
	}
	return P / (V * 4);
}

/*inline dmat3 full_mat3(double a) {
    const double m[] = { a, a, a, a, a, a, a, a, a};
    return glm::make_mat3(m);
}

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
dmat3 moment_of_inertia(const imesh3& mesh) {
	constexpr double a = 1 / 60., b = 1 / 120.;
    constexpr double canon[] = { a, b, b, b, a, b, b, b, a};
	const dmat3 canonical = glm::make_mat3(canon);
    dmat3 C = dmat3(0); // covariance
	for (const auto& f : mesh) {
		dmat3 A;
		for (auto i : range(3))
			glm::row(A, i) = f[i]; // TODO: or column(i) = ?
		C += (transpose(A) * canonical * A) * determinant(A);
	}
	return full_mat3(C[0][0] + C[1][1] + C[2][2]) - C; // C -> I
}*/

bool is_aabb(const imesh3& mesh) {
	aabb<ivec3> box(mesh);

	// All vertices must be made from extreme coordinates
	for (auto f : mesh)
		for (ivec3 v : f)
			if (v.x != box.min[0] && v.x != box.max[0]
			 && v.y != box.min[1] && v.y != box.max[1]
			 && v.z != box.min[2] && v.z != box.max[2])
				return false;

	// Every face must have one coordinate constant
	for (auto f : mesh) {
		bool xx = f.a.x == f.b.x && f.b.x == f.c.x;
		bool yy = f.a.y == f.b.y && f.b.y == f.c.y;
		bool zz = f.a.z == f.b.z && f.b.z == f.c.z;
		if (!xx && !yy && !zz)
			return false;
	}
	return true;
}
