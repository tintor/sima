#include "properties.h"
#include "aabb.h"
#include "scalar.h"
#include "primitives.h"

double triangle_volume(double3 a, double3 b, double3 c) {
	double z = a.z;
	z += b.z;
	z += c.z;
	double s;
	s = (a.y + b.y) * (a.x - b.x);
	s += (b.y + c.y) * (b.x - c.x);
	s += (c.y + a.y) * (c.x - a.x);
	return s * z;
}

double signed_volume(double3 a, double3 b, double3 c, double3 d) {
	double v;
	v = triangle_volume(a, b, c);
	v += triangle_volume(d, b, a);
	v += triangle_volume(a, c, b);
	v += triangle_volume(d, a, c);
	return v / 6.0;
}

double signed_volume(const mesh3& mesh) {
	double v = 0;
	for (const triangle3& f : mesh)
		v += triangle_volume(f.a, f.b, f.c);
	return v / 6.0;
}

double volume(const mesh3& mesh) {
	return std::abs(signed_volume(mesh));
}

// TODO test with volume of cube (randomly rotated)

// Center of mass of a valid polyhedron
double3 center_of_mass(const mesh3& mesh) {
	double3 P = {0, 0, 0};
	double V = 0;
	for (const triangle3& f : mesh) {
		double v = dot(f.a, cross(f.b, f.c));
		P += (f.a + f.b + f.c) * v;
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

bool is_aabb(const mesh3& mesh) {
	aabb3 box(mesh);

	// All vertices must be made from extreme coordinates
	for (auto f : mesh)
		for (double3 v : f)
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
