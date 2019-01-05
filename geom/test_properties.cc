#include <geom/properties.h>
#include <geom/generators.h>
#include <catch.hpp>
#include <core/align_alloc.h>

// TODO center of mass
// TODO translate mesh
// TODO rotate mesh
TEST_CASE("volume_of_cube", "[properties]") {
	auto m = generate_box(2, 3, 5);
	REQUIRE(SignedVolume(m) == 240);
}

double min_eigenvalue(double xx, double yy, double zz, double xy, double xz, double yz) {
	double a = xx + yy + zz;
	double b = xx*zz + yy*zz + xx*yy - yz*yz - xy*xy - xz*xz;
	double c = xx*yy*zz + 2*xy*xz*yz - xz*xz*yy - xy*xy*zz - xx*yz*yz;
	a /= 3;
	REQUIRE(a*a >= b/3);
	double q = a*a - b/3;
	double q_sqrt = sqrt(q);
	double t = (a*(b/2 - a*a) - c/2) / (q_sqrt * q);
	REQUIRE(t <= 1.0);
	REQUIRE(t >= -1.0);
    return a - 2 * q_sqrt * cos(acos(t) / 3);
}

/*	auto ee = min_eigenvalue(ss.x, ss.y, ss.z, xy, xz, yz);
	print("act %s\n", ee);*/

/*	auto ev = eigenvector(ss.x - ee, ss.y - ee, ss.z - ee, xy, xz, yz);
	print("(%s)\n", ev);
	// M=(A-I*e), where e is minimal eigen value
	// solve v given M*v=0*/

// Test cases:
// - all identical points
// - all points on one line
// - all points uniform in cube
// - all points on sphere surface
TEST_CASE("eigen - plane normal", "[eigen]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	aligned_vector<double3> points;
	points.resize(20);
	for (int j = 0; j < 200000; j++) {
		double3 normal = uniform_dir3(rnd);
		for (int i = 0; i < points.size(); i++) {
			double3 a = uniform3(rnd, -1, 1);
			double aa = dot(a, normal);
			a -= aa * normal;
			double eps = 1e-4;
			a += uniform(rnd, -eps, +eps) * normal;
			points[i] = a;
		}
		double3 n = eigen_vector(points);
		double dn = dot(normal, n);
		REQUIRE(abs(dn) >= 1 - 1e-6);
	}
}
