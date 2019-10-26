#include <sim/gauss_seidel.h>
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
