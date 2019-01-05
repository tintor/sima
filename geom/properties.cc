#include <Eigen/Dense>
#include <geom/properties.h>
#include <geom/aabb.h>
#include <geom/matrix.h>
#include <geom/primitives.h>

double SignedTriangleVolume6(double3 a, double3 b, double3 c) {
	double z = a.z;
	z += b.z;
	z += c.z;
	double s;
	s = (a.y + b.y) * (a.x - b.x);
	s += (b.y + c.y) * (b.x - c.x);
	s += (c.y + a.y) * (c.x - a.x);
	return s * z;
}

double SignedRingVolume6(cspan<double3> ring) {
	double z = 0;
	for (auto e : ring)
		 z += e.z;
	double s = 0;
	for (auto e : Edges(ring))
		s += (e.a.y + e.b.y) * (e.a.x - e.b.x);
	return s * z;
}

double SignedVolume(const mesh3& mesh) {
	double v = 0;
	for (const triangle3& f : mesh)
		v += SignedTriangleVolume6(f.a, f.b, f.c);
	return v / 6.0;
}

double SignedVolume(const xmesh3& mesh) {
	double v = 0;
	for (face f : mesh)
		for (auto ring : f)
		v += SignedRingVolume6(ring);
	return v / 6.0;
}

// TODO test with volume of cube (randomly rotated)

// Center of mass of a valid polyhedron
double3 CenterOfMass(const mesh3& mesh) {
	double3 P = {0, 0, 0};
	double V = 0;
	for (const triangle3& f : mesh) {
		double v = dot(f.a, cross(f.b, f.c));
		P += (f.a + f.b + f.c) * v;
		V += v;
	}
	return P / (V * 4);
}

double2 centroid(cspan<double2> poly) {
	double2 P = {0, 0};
	double2 V = 0;
	for (auto i : range(poly.size())) {
		double2 a = poly[i];
		double2 b = poly[(i + 1) % poly.size()];
		double v = a.x * b.y - a.y * b.x;
		P += (a + b) * v;
		V += v;
	}
	return P / (V * 3);
}

double2 CenterOfMass(const polygon2& poly) {
	double2 P = {0, 0};
	double2 V = 0;
	for (auto [a, b] : Edges(poly)) {
		double v = a.x * b.y - a.y * b.x;
		P += (a + b) * v;
		V += v;
	}
	return P / (V * 3);
}

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
double33 moment_of_inertia(const mesh3& mesh) {
	constexpr double a = 1 / 60., b = 1 / 120.;
	const double33 canonical = {{a, b, b}, {b, a, b}, {b, b, a}};
    auto C = double33::make(0); // covariance
	for (const auto& f : mesh) {
		double33 A;
		// TODO setting rows or columns?
		A.a = f[0];
		A.b = f[1];
		A.c = f[2];
		C += transpose(A) * canonical * A * det(A);
	}
	return double33::make(trace(C)) - C; // C -> I
}

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

double3 eigen_vector(cspan<double3> points) {
	// compute mean
	double3 m = {0, 0, 0};
	for (auto p : points)
		m += p;
	m /= points.size();

	// compute matrix
	double3 ss = {0, 0, 0};
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
	return {vectors(0, i).real(), vectors(1, i).real(), vectors(2, i).real()};
}
