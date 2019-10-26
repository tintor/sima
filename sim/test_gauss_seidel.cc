#include <sim/gauss_seidel.h>
#include <sim/newton.h>
#include <catch.hpp>

TEST_CASE("gauss seidel basic", "[gauss_seidel]") {
	vector<double> a = { 16, 3, 7, -11 };
	vector<double> b = { 11, 13 };
	vector<double> x = { 0, 0 };
	vector<double> y = { 0, 0 };

	print("b = %.10f\n", b);
	for (auto i : range(20)) {
		GaussSeidelSolver(a, b, 1, x);
		y[0] = dot(cspan<double>(a).subspan(0, 2), x);
		y[1] = dot(cspan<double>(a).subspan(2, 2), x);
		print("%s: x = %.10f, A * x = %.10f\n", i, x, y);
	}

	REQUIRE(b[0] == Approx(y[0]).margin(1e-10));
	REQUIRE(b[1] == Approx(y[1]).margin(1e-10));
}

TEST_CASE("newton sqrt()", "[newton]") {
	for (double e : { 2, 3, 4, 5, 6, 7, 8, 9, 10 }) {
		double a = NewtonMethod(1, 5, [e](double x) { return x * x - e; }, [](double x) { return 2 * x; });
		REQUIRE(abs(sqrt(e) - a) < 1e-8);
	}
}
