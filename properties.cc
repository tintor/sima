#include "properties.h"
#include "aabb.h"

// Volume of a valid polyhedron
double volume(const imesh3& mesh) {
	double volume = 0;
	for (const auto& f : mesh) {
		long s = 0, z = 0;
		for (auto [a, b] : f.edges()) {
			z += b.z;
			s += (long(a.y) + b.y) * (long(a.x) - b.x);
		}
		volume += double(z) * s;
	}
	return volume / 6;
}

// TODO test with volume of cube (randomly rotated)

// Center of mass of a valid polyhedron
ivec3 center_of_mass(const imesh3& mesh) {
	dvec3 P = {0, 0, 0};
	double V = 0;
	for (const auto& f : mesh) {
		double v = dot(f.a, cross(f.b, f.c));
		P += dvec3(f.a + f.b + f.c) * v;
		V += v;
	}
	P /= V * 4;
	return {std::round(P.x), std::round(P.y), std::round(P.z)};
}

inline dmat3 full_mat3(double a) {
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
}

bool is_aabb(const imesh3& mesh) {
	aabb<ivec3> box(mesh);

	// All vertices must be made from extreme coordinates
	for (auto f : mesh)
		for (ivec3 v : f)
			if (v.x != box.mm[0].min && v.x != box.mm[0].max
			 && v.y != box.mm[1].min && v.y != box.mm[1].max
			 && v.z != box.mm[2].min && v.z != box.mm[2].max)
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
