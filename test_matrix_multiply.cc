#include <core/format.h>
#include <core/range.h>
#include <core/timestamp.h>
#include <immintrin.h>

#include <catch.hpp>

template <int N>
struct Matrix {
    std::array<float, N * N> m;

    const float& operator()(uint r, uint c) const { return m[r * N + c]; }
    float& operator()(uint r, uint c) { return m[r * N + c]; }
};

template <int N>
void multiply(const Matrix<N>& a, const Matrix<N>& b, Matrix<N>& c) {
    for (uint row = 0; row < N; row++) {
        for (uint col = 0; col < N; col++) c(row, col) = 0;
        for (uint i = 0; i < N; i++)
            for (uint col = 0; col < N; col++) c(row, col) += a(row, i) * b(i, col);
    }
}

inline uint min(uint a, uint b) { return (a < b) ? a : b; }

template <int N>
void multiply2(const Matrix<N>& a, const Matrix<N>& b, Matrix<N>& c) {
    const uint s = 512;  // depends on cache size

    memset(c.m.data(), 0, 4 * N * N);

    for (uint ki = 0; ki < N; ki += s)
        for (uint kc = 0; kc < N; kc += s)
            for (uint row = 0; row < N; row++)
                for (uint i = ki; i < min(ki + s, N); i++)
                    for (uint col = kc; col < min(kc + s, N); col++) c(row, col) += a(row, i) * b(i, col);
}

TEST_CASE("fma flops/cycle", "[!hide][fma]") {
    std::default_random_engine rnd;
    float8 a, b, c;
    for (int i = 0; i < 8; i++) a[i] = std::uniform_real_distribution(-1., 1.)(rnd);
    for (int i = 0; i < 8; i++) b[i] = std::uniform_real_distribution(-1., 1.)(rnd);
    for (int i = 0; i < 8; i++) c[i] = std::uniform_real_distribution(-1., 1.)(rnd);

    float af = a.x;
    float bf = b.x;
    float cf = c.x;

    const ulong N = 10000000000;
    Timestamp t0;
    for (auto i : range(N)) {
        a = _mm256_fmadd_ps(a, b, c);
    }
    Timestamp t1;
    float s = 0;
    for (int i = 0; i < 8; i++) s += a[i];
    print("%f %f %s\n", s, a.x, (double(N) * 16) / t0.elapsed(t1));

    Timestamp t2;
    for (auto i : range(N)) {
        af = af * bf + cf;
    }
    Timestamp t3;
    print("%f %s\n", af, (double(N) * 2) / t2.elapsed(t3));
}

TEST_CASE("matrix multitply benchmark", "[!hide][matrix_multiply]") {
    std::default_random_engine rnd;
    auto* a = new Matrix<2048>;
    auto* b = new Matrix<2048>;
    auto* c = new Matrix<2048>;

    for (float& e : c->m) e = std::uniform_real_distribution(-1., 1.)(rnd);
    for (float& e : a->m) e = std::uniform_real_distribution(-1., 1.)(rnd);
    for (float& e : b->m) e = std::uniform_real_distribution(-1., 1.)(rnd);
    const int N = 2048;

    SECTION("kernel") {
        Timestamp t0;
        multiply(*a, *b, *c);
        Timestamp t1;
        print("flops/cycle %s\n", (double(N) * double(N) * N * 2) / t0.elapsed(t1));
    }
    SECTION("kernel2") {
        Timestamp t0;
        multiply2(*a, *b, *c);
        Timestamp t1;
        print("flops/cycle %s\n", (double(N) * double(N) * N * 2) / t0.elapsed(t1));
    }

    float s = 0;
    for (float& e : c->m) s += e;
    print("%f\n", s);
}
