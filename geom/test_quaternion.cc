#include <catch.hpp>
#include <geom/quaternion.h>
#include <core/range.h>

TEST_CASE("quat from axis angle", "[quaternion]") {
	quat q = quat_from_axis_angle(double4{1, 0, 0}, PI / 2);
	auto v = quat_rotate(q, double3{0, 1, 0});
	REQUIRE(length(v - double3{0, 0, 1}) <= 1e-15);
}

TEST_CASE("quat rotate", "[quaternion]") {
	std::default_random_engine rnd;
	for (auto i : range(10000)) {
		quat q = uniform_dir4(rnd);
		double3 a = uniform3(rnd, 0, 1).xyz;

		double3 b = quat_rotate(q, a);
		double3 c = mul(quat_to_matrix(q), a);
		REQUIRE(length(b - c) < 1e-15);
	}
}

vector<pair<quat, double3>> g_quat_vec;
vector<double3> g_out_vec;

class Init {
public:
	Init() {
		std::default_random_engine rnd;
		constexpr int N = 1000000;
		g_quat_vec.resize(N);
		g_out_vec.resize(N);
		for (auto i : range(N))
			g_quat_vec[i] = { uniform_dir4(rnd), uniform3(rnd, 0, 1).xyz };
	}
} g_init;

TEST_CASE("quat rotate perf", "[quaternion][!hide]") {
	const auto n = g_quat_vec.size();
	for (auto k : range(100)) {
		for (auto i : range(n)) {
			auto [q, v] = g_quat_vec[i];
			g_out_vec[i] = quat_rotate(q, v);
		}
	}
}

TEST_CASE("quat matrix perf", "[quaternion][!hide]") {
	auto m = quat_to_matrix(g_quat_vec[0].first);
	const auto n = g_quat_vec.size();
	for (auto k : range(100)) {
		for (auto i : range(n)) {
			auto v = g_quat_vec[i].second;
			g_out_vec[i] = mul(m, v);
		}
	}
}
