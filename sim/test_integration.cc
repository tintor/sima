#include <core/range.h>
#include <sim/integration.h>

#include <catch.hpp>

TEST_CASE("Numerical integrators") {
    // differential edquation: y' = y, y(0) = 1
    // exact solution: y(t) = e^t
    // simulating 320 steps from 0.0 to 4.0, step size 0.0125
    const auto& der = [](double e, double t) { return e; };
    double k = 4.0 / 320;

    double e = 1;
    for (auto i : range(320)) e = EulersMethod<double, double>(e, i * k, k, der);
    CHECK(abs(e - exp(4.0)) <= 1.34);

    e = 1;
    for (auto i : range(320)) e = MidpointMethod<double, double>(e, i * k, k, der);
    CHECK(abs(e - exp(4.0)) <= 0.0057);

    e = 1;
    for (auto i : range(320)) e = RungeKutta4<double, double>(e, i * k, k, der);
    CHECK(abs(e - exp(4.0)) <= 4.4e-08);
}
