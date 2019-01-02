#include <catch.hpp>
#include <core/exception.h>
#include <core/range.h>
#include <geom/matrix.h>

TEST_CASE("inv(double22)", "[matrix]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	for (auto i : range(100)) {
		double22 m;
		m.a = uniform2(rnd, 0, 1);
		m.b = uniform2(rnd, 0, 1);

		double22 n = inv(m);
		double22 e = mul(m, n);
		REQUIRE_NEAR(e.a, d2(1, 0), 5e-14);
		REQUIRE_NEAR(e.b, d2(0, 1), 5e-14);
	}
}

TEST_CASE("inv(double33)", "[matrix]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	for (auto i : range(100)) {
		double33 m;
		m.a = uniform3(rnd, 0, 1);
		m.b = uniform3(rnd, 0, 1);
		m.c = uniform3(rnd, 0, 1);

		double33 n = inv(m);
		double33 e = mul(m, n);
		REQUIRE_NEAR(e.a, d3(1, 0, 0), 1e-13);
		REQUIRE_NEAR(e.b, d3(0, 1, 0), 1e-13);
		REQUIRE_NEAR(e.c, d3(0, 0, 1), 1e-13);
	}
}

TEST_CASE("inv(double44)", "[matrix]") {
	std::default_random_engine rnd;
	rnd.seed(0);
	for (auto i : range(100)) {
		double44 m;
		m.a = uniform4(rnd, 0, 1);
		m.b = uniform4(rnd, 0, 1);
		m.c = uniform4(rnd, 0, 1);
		m.d = uniform4(rnd, 0, 1);

		double44 n = inv(m);
		double44 e = mul(m, n);
		REQUIRE_NEAR(e.a, d4(1, 0, 0, 0), 3e-13);
		REQUIRE_NEAR(e.b, d4(0, 1, 0, 0), 3e-13);
		REQUIRE_NEAR(e.c, d4(0, 0, 1, 0), 3e-13);
		REQUIRE_NEAR(e.d, d4(0, 0, 0, 1), 3e-13);
	}
}
