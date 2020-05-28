#pragma once

#include <core/range.h>
#include <core/std.h>

inline double dot(cspan<double> a, cspan<double> b) {
    // TODO simd, do 4 loop iterations at once
    double sum = 0;
    for (auto i : range(a.size())) sum += a[i] * b[i];
    return sum;
}

// Solver for system of linear equations:
// A - n x n matrix
// b - n x 1 vector
// x - n x 1 vector
// Find x such that A * x = b
inline void GaussSeidelSolver(cspan<double> a, cspan<double> b, int iterations, vector<double>& x) {
    auto n = b.size();
    assert(a.size() == n * n);
    assert(x.size() == n);

    for (auto it : range(iterations))
        for (auto i : range(n)) {
            auto ai = a.subspan(i * n, n);
            double sigma = dot(ai, x) - ai[i] * x[i];
            x[i] = (b[i] - sigma) / ai[i];
        }
}

// Successive over-relaxation variant
inline void GaussSeidelSolverSOR(cspan<double> a, cspan<double> b, double w, int iterations, vector<double>& x) {
    auto n = b.size();
    assert(a.size() == n * n);
    assert(x.size() == n);

    for (auto it : range(iterations))
        for (auto i : range(n)) {
            auto ai = a.subspan(i * n, n);
            double sigma = dot(ai, x) - ai[i] * x[i];
            x[i] += w * ((b[i] - sigma) / ai[i] - x[i]);
        }
}
