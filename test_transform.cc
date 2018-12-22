#include "transform.h"
#include <core/timestamp.h>
#include "catch.hpp"

/*TEST_CASE("translate", "[transform]") {
	ivec4 t = {1, 2, 3, 0};
	vector<ivec4> a;
	a.resize(1024 + 1);

	for (uint i = 0; i < 1024; i++) {
		// init a
		for (uint j = 0; j < i + 1; j++)
			a[j] = {10, 20, 30, 0};

		translate(span<ivec4>(a.data(), a.data() + i), t);

		// verify a
		for (uint j = 0; j < i; j++) {
			REQUIRE(a[j].x == 11);
			REQUIRE(a[j].y == 22);
			REQUIRE(a[j].z == 33);
			REQUIRE(a[j].w == 0);
		}
		REQUIRE(a[i].x == 10);
		REQUIRE(a[i].y == 20);
		REQUIRE(a[i].z == 30);
		REQUIRE(a[i].w == 0);
	}
}

Timestamp g_ts;

void reset(string_view name = "") {
	if (name.size() > 0)
		print("%s: %.3lf\n", name, g_ts.elapsed_s());
	g_ts = Timestamp();
}

TEST_CASE("translate_simd benchmark", "[.][transform][benchmark]") {
	ivec4 t = {1, 2, 3, 0};
	vector<ivec4> a;
	a.resize(1 << 24); // 256MB

	for (ivec4& v : a)
		v = {10, 20, 30, 0};

	Timestamp ta;
	translate(a, t);
	Timestamp tb;

	double seconds = ta.elapsed_s(tb);
	print("%s ivec4s translated in %.3f seconds (%.0f vectors per second)",
		a.size(), seconds, a.size() / seconds);
}*/
