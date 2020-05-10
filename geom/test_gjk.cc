#include <core/exception.h>
#include <geom/gjk.h>
#include <geom/vector.h>

#include <catch.hpp>

template <typename RNG, typename Convex>
void verify_support(RNG& rnd, const Convex& support) {
    for (auto i : range(1000000)) {
        double2 d1 = uniform_dir2(rnd);
        double2 d2 = uniform_dir2(rnd);
        double2 s1 = support(d1);
        double2 s2 = support(d2);
        REQUIRE(dot(d1, s1) >= dot(d1, s2));
        REQUIRE(dot(d2, s2) >= dot(d2, s1));
    }
}

template <typename ConvexA, typename ConvexB>
bool brute2_intersects(const ConvexA& a, const ConvexB& b, uint iterations = 1000000) {
    std::default_random_engine rnd;
    rnd.seed(0);
    for (auto i : range(iterations)) {
        double2 dir = uniform_dir2(rnd);
        if (i == 0) dir = {-1, -1};
        if (i == 1) dir = {1, 1};
        auto sa = a(dir);
        auto sb = b(-dir);
        auto s = sa - sb;
        if (dot(s, dir) < 0) return false;
    }
    return true;
}

template <typename RNG, typename ConvexA, typename ConvexB>
optional<pair<double3, double3>> brute3(RNG& rnd, const ConvexA& a, const ConvexB& b, uint iterations = 1000000) {
    optional<pair<double3, double3>> res;
    double m = -INF;
    for (auto i : range(iterations)) {
        double3 dir = uniform_dir3(rnd);
        auto sa = a(dir);
        auto sb = b(-dir);
        auto s = sa - sb;
        if (dot(s, dir) < 0) {
            double cm = length(s);
            if (cm > m) {
                m = cm;
                res = {sa, sb};
            }
        }
    }
    return (m >= 0) ? res : nullopt;
}

// TODO test correctness for curved shapes and linear shapes
// TODO verify number of iterations for linear shapes (there should be never more than dim(V) + 1 iterations)
// TODO verify with different random dirs
// TODO verify there will be only one iteration if prev dir is used

bool eq(pair<double2, double2> a, pair<double2, double2> b) {
    return all(a.first == b.first) && all(a.second == b.second);
}

bool eq(optional<pair<double2, double2>> a, optional<pair<double2, double2>> b) {
    return (a.has_value() && b.has_value() && eq(*a, *b)) || (!a.has_value() && !b.has_value());
}

TEST_CASE("gjk2 point/triangle", "[.][gjk2]") {
    double2 a = {1, 1};
    double2 b = {5, 1};
    double2 c = {1, 4};
    auto out = [=](double2 p) {
        int counter = 0;
        auto support_a = [a, b, c, &counter](double2 d) {
            counter += 1;
            double aa = dot(a, d);
            double bb = dot(b, d);
            double cc = dot(c, d);
            double m = max(aa, bb, cc);
            REQUIRE(m >= aa);
            REQUIRE(m >= bb);
            REQUIRE(m >= cc);
            if (m == aa) return a;
            if (m == bb) return b;
            if (m == cc) return c;
            THROW(runtime_error);
        };
        auto support_b = [p](double2 d) { return p; };

        double2 dir = {1, 1};
        auto res = gjk_classify(support_a, support_b, dir) >= 0;
        int snapshot = counter;
        print("gjk([%s] [%s] [%s] vs [%s]) -> %s (counter %s)\n", a, b, c, p, res, snapshot);

        auto res2 = brute2_intersects(support_a, support_b);
        print("brute([%s] [%s] [%s] vs [%s]) -> %s\n\n", a, b, c, p, res2);

        REQUIRE(res == res2);
    };
    out(double2{0, 0});
    out(double2{1, 0});
    out(double2{0, 1});
    out(double2{9, -2});
    out(a);
    out(b);
    out(c);
    out((a + b) / 2);
    out((a + c) / 2);
    out((b + c) / 2);
}

TEST_CASE("gjk2 circle/circle", "[.][gjk]") {
    std::default_random_engine rnd;
    rnd.seed(0);
    int m = 0;
    for (int j : range(20)) {
        for (int i : range(-10, 10)) {
            if (i == 0) continue;
            double delta = sign(i) / pow(10, abs(i));
            double2 a = {5, 1};
            double2 b = {1, 4};
            int counter = 0;
            auto support_a = [a, delta, &counter](double2 d) {
                counter += 1;
                return normalize(d) * (3 + delta) + a;
            };
            auto support_b = [b, delta](double2 d) { return normalize(d) * (2 + delta) + b; };
            double2 dir = uniform_dir2(rnd);
            bool actual = gjk_classify(support_a, support_b, dir, 1000) >= 0;
            maximize(m, counter);
            bool expected = brute2_intersects(support_a, support_b);
            ASSERT_ALWAYS(actual == expected, "i=%d", i);
        }
    }
    print("max counter %s\n", m);
}

// generate random triangle / tetrahedron cases
// generate random sphere / sphere cases
// generate random capsule / capsule cases
// cross verify with (polygon / polyhedron) (distance / intersects) functions
