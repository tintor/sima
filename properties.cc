#include "properties.h"
#include "aabb.h"
#include "scalar.h"
#include "ivec.h"

int128 signed_volume_mul6(const imesh3& mesh) {
	int128 volume = 0;
	for (const itriangle3& f : mesh) {
		int128 s = 0;
		long z = 0;
		for (auto [a, b] : f.edges()) {
			z = addi(z, b.z);
			long y = addi(a.y, b.y);
			long x = subi(a.x, b.x);
			s = addi(s, muli(y, x));
		}
		int128 v = muli(s, z);
		volume = addi(volume, v);
	}
	return volume;
}

double volume(const imesh3& mesh) {
	return std::abs(signed_volume_mul6(mesh)) / 6.0;
}

// TODO test with volume of cube (randomly rotated)

// Center of mass of a valid polyhedron
ivec3 center_of_mass(const imesh3& mesh) {
	lvec3 P = {0, 0, 0};
	long V = 0;
	for (const auto& f : mesh) {
		long v = doti(f.a, crossi(f.b, f.c));
		lvec3 p = muli(addi(addi(f.a, f.b), f.c), v);
		P = addi(P, p);
		V = addi(V, v);
	}
	return vconvert(P / muli(V, 4), ivec3);
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
