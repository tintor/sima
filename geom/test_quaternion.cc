#include <core/exception.h>
#include <core/range.h>
#include <geom/quaternion.h>

#include <catch.hpp>

TEST_CASE("identity quat", "[quaternion]") {
    double3 v = {1, 2, 3};
    double3 r = quat_rotate(IDENTITY_QUAT, v);
    REQUIRE_NEAR(v, r, 1e-15);

    double33 m = quat_to_matrix(IDENTITY_QUAT);
    REQUIRE_NEAR(m.a, d3(1, 0, 0), 1e-15);
    REQUIRE_NEAR(m.b, d3(0, 1, 0), 1e-15);
    REQUIRE_NEAR(m.c, d3(0, 0, 1), 1e-15);
}

TEST_CASE("identity mat", "[quaternion]") {
    double3 v = {1, 2, 3};
    double33 m = {double3{1, 0, 0}, double3{0, 1, 0}, double3{0, 0, 1}};
    REQUIRE_NEAR(v, m * v, 1e-15);
}

TEST_CASE("mat mul", "[quaternion]") {
    double a = PI / 3;
    double c = cos(a), s = sin(a);
    double33 m = {double3{c, s, 0}, double3{-s, c, 0}, double3{0, 0, 1}};

    double3 v = {2, 0, 0};
    double3 w = {2 * cos(PI / 3), 2 * sin(PI / 3), 0};
    REQUIRE_NEAR(w, m * v, 1e-15);
}

TEST_CASE("quat from axis angle", "[quaternion]") {
    quat q = quat_from_axis_angle(double3{1, 0, 0}, PI / 2);
    auto v = quat_rotate(q, double3{0, 1, 0});
    REQUIRE_NEAR(v, d3(0, 0, 1), 1e-15);
}

TEST_CASE("quat and matrix", "[quaternion]") {
    std::default_random_engine rnd;
    rnd.seed(0);
    for (auto i : range(10000)) {
        quat q = uniform_dir4(rnd);
        REQUIRE_NEAR(length(q), 1, 1e-15);

        double33 m = quat_to_matrix(q);
        double33 mm = transpose(quat_to_matrix(quat_inv(q)));
        REQUIRE_NEAR(m.a, mm.a, 1e-15);
        REQUIRE_NEAR(m.b, mm.b, 1e-15);
        REQUIRE_NEAR(m.c, mm.c, 1e-15);

        REQUIRE_NEAR(dot(m.a, m.b), 0, 1e-15);
        REQUIRE_NEAR(dot(m.a, m.c), 0, 1e-15);
        REQUIRE_NEAR(dot(m.b, m.c), 0, 1e-15);

        REQUIRE_NEAR(length(m.a), 1, 1e-14);
        REQUIRE_NEAR(length(m.b), 1, 1e-14);
        REQUIRE_NEAR(length(m.c), 1, 1e-14);

        REQUIRE_NEAR(det(m), 1, 1e-14);

        // TODO error is huge!
        quat qq = quat_from_matrix(m);
        REQUIRE((length(q - qq) <= 1e-11 || length(q + qq) <= 1e-11));
    }
}

TEST_CASE("quat slerp", "[quaternion]") {
    std::default_random_engine rnd;
    rnd.seed(0);
    for (auto i : range(10000)) {
        quat a = uniform_dir4(rnd);
        double step = M_PI / 5;
        quat s = quat_from_axis_angle(uniform_dir3(rnd), step);
        quat a1 = quat_mul(a, s);
        quat a2 = quat_mul(a1, s);
        quat a3 = quat_mul(a2, s);
        quat a4 = quat_mul(a3, s);
        REQUIRE_NEAR(a1, slerp(a, a2, 0.5), 1e-15);
        REQUIRE_NEAR(a1, slerp(a, a3, 1. / 3), 1e-15);
        REQUIRE_NEAR(a2, slerp(a, a3, 2. / 3), 1e-15);
        REQUIRE_NEAR(a1, slerp(a, a4, 0.25), 1e-15);
    }
}

TEST_CASE("quat rotate", "[quaternion]") {
    std::default_random_engine rnd;
    rnd.seed(0);
    for (auto i : range(10000)) {
        quat q = uniform_dir4(rnd);
        double3 a = uniform3(rnd, 0, 1).xyz;

        double3 b = quat_rotate(q, a);
        double3 c = quat_to_matrix(q) * a;
        REQUIRE_NEAR(b, c, 1e-15);
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
        for (auto i : range(N)) g_quat_vec[i] = {uniform_dir4(rnd), uniform3(rnd, 0, 1).xyz};
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
            g_out_vec[i] = m * v;
        }
    }
}
