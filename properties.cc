#include "properties.h"
#include "aabb.h"
#include "scalar.h"
#include "primitives.h"
#include <Eigen/Dense>

double triangle_volume(double4 a, double4 b, double4 c) {
	double z = a.z;
	z += b.z;
	z += c.z;
	double s;
	s = (a.y + b.y) * (a.x - b.x);
	s += (b.y + c.y) * (b.x - c.x);
	s += (c.y + a.y) * (c.x - a.x);
	return s * z;
}

double signed_volume(double4 a, double4 b, double4 c, double4 d) {
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
double4 center_of_mass(const mesh3& mesh) {
	double4 P = {0, 0, 0};
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
	aabb4 box(mesh);

	// All vertices must be made from extreme coordinates
	for (auto f : mesh)
		for (double4 v : f)
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

double4 eigen_vector(span<const point3> points) {
	// compute mean
	double4 m = {0, 0, 0, 0};
	for (point3 p : points)
		m += p;
	m /= points.size();

	// compute matrix
	double4 ss = {0, 0, 0, 0};
	double xy = 0, xz = 0, yz = 0;
	for (auto p : points) {
		p -= m;
		ss += p * p;
		xy += p.x * p.y;
		xz += p.x * p.z;
		yz += p.y * p.z;
	}

	Eigen::MatrixXd a = Eigen::MatrixXd::Zero(3, 3);
	a(0, 0) = ss.x;
	a(1, 1) = ss.y;
	a(2, 2) = ss.z;
	a(0, 1) = a(1, 0) = xy;
	a(0, 2) = a(2, 0) = xz;
	a(1, 2) = a(2, 1) = yz;

	Eigen::EigenSolver<Eigen::MatrixXd> es(a, true);
	auto values = es.eigenvalues();
	int i = 0;
	if (values(1).real() < values(i).real())
		i = 1;
	if (values(2).real() < values(i).real())
		i = 2;

	auto vectors = es.eigenvectors();
	return {vectors(0, i).real(), vectors(1, i).real(), vectors(2, i).real(), 0};
}
