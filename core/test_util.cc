#include <core/align_alloc.h>
#include <core/format.h>
#include <core/range.h>
#include <core/string_util.h>
#include <core/util.h>
#include <immintrin.h>

#include <catch.hpp>
#include <range/v3/algorithm.hpp>

TEST_CASE("ranges basic", "[ranges]") {
    vector<int> v = {6, 2, 3, 4, 5, 6};
    cout << "vector:   " << ranges::count_if(v, [](int i) { return i == 6; }) << endl;
}

TEST_CASE("min", "[util]") {
    REQUIRE(min(2, 1) == 1);
    REQUIRE(min(5, 7, 9) == 5);
    REQUIRE(min(9, 7, 5) == 5);
    REQUIRE(min(3, 0, -1, 6, 10, -4, 3) == -4);
}

// Two changes to reduce numerical error:
// - double as accumulator even if T is float
// - recursive subdivision implementation reduces numerical error to O(log(N)) from O(N) of naive algorithm
template <typename T>
double FSum(cspan<T> v) {
    if (v.size() >= 4) return FSum(v.first(v.size() / 2)) + FSum(v.pop_front(v.size() / 2));
    if (v.size() == 3) return v[0] + v[1] + v[2];
    if (v.size() == 2) return v[0] + v[1];
    if (v.size() == 1) return v[0];
    return 0;
}

template <typename T>
double M2(cspan<T> v, double mean) {
    if (v.size() >= 4) return M2(v.first(v.size() / 2), mean) + M2(v.pop_first(v.size() / 2), mean);
    if (v.size() == 3) return sqr(v[0] - mean) + sqr(v[1] - mean) + sqr(v[2] - mean);
    if (v.size() == 2) return sqr(v[0] - mean) + sqr(v[1] - mean);
    if (v.size() == 1) return sqr(v[0] - mean);
    return 0;
}

template <typename T>
double Mean(cspan<T> v) {
    return FSum(v) / v.size();
}

template <typename T>
double Variance(cspan<T> v, double mean) {
    return M2(v, mean) / v.size();
}

template <typename T>
T NaiveSum(cspan<float> v) {
    T a = 0;
    for (auto e : v) a += e;
    return a;
}

template <typename T>
T NaiveSum(cspan<double> v) {
    T a = 0;
    for (auto e : v) a += e;
    return a;
}

template <typename T>
T ShardedSum2(cspan<float> v) {
    auto n = v.size() / 8;
    return NaiveSum<T>(v.subspan(n * 0, n)) + NaiveSum<T>(v.subspan(n * 1, n)) + NaiveSum<T>(v.subspan(n * 2, n)) +
           NaiveSum<T>(v.subspan(n * 3, n)) + NaiveSum<T>(v.subspan(n * 4, n)) + NaiveSum<T>(v.subspan(n * 5, n)) +
           NaiveSum<T>(v.subspan(n * 6, n)) + NaiveSum<T>(v.subspan(n * 7, n));
}

template <typename T>
T ShardedSum2(cspan<double> v) {
    auto n = v.size() / 8;
    return NaiveSum<T>(v.subspan(n * 0, n)) + NaiveSum<T>(v.subspan(n * 1, n)) + NaiveSum<T>(v.subspan(n * 2, n)) +
           NaiveSum<T>(v.subspan(n * 3, n)) + NaiveSum<T>(v.subspan(n * 4, n)) + NaiveSum<T>(v.subspan(n * 5, n)) +
           NaiveSum<T>(v.subspan(n * 6, n)) + NaiveSum<T>(v.subspan(n * 7, n));
}

float ShardedSumVectorizedF(cspan<float> v) {
    float4 a = {0, 0, 0, 0};
    for (size_t i = 0; i < v.size(); i += 4) a += *(float*)(v.data() + i);
    return a[0] + a[1] + a[2] + a[3];
}

float ShardedSumSimdF(cspan<float> v) {
    float4 a = {0, 0, 0, 0};
    for (size_t i = 0; i < v.size(); i += 4) a += _mm_load_ps(v.data() + i);
    return a[0] + a[1] + a[2] + a[3];
}

template <typename T>
T ShardedSum(cspan<float> v) {
    T a0 = 0;
    T a1 = 0;
    T a2 = 0;
    T a3 = 0;
    for (size_t i = 0; i < v.size(); i += 4) {
        a0 += v[i];
        a1 += v[i + 1];
        a2 += v[i + 2];
        a3 += v[i + 3];
    }
    return a0 + a1 + a2 + a3;
}

template <typename T>
T IncrementalMean(cspan<float> v) {
    size_t count = 0;
    T mean = 0;
    for (float a : v) mean += (a - mean) / ++count;
    return mean;
}

#define L2(A) function<double()>([&]() { return double(A); })

TEST_CASE("sum", "[util][.]") {
    for (size_t seed = 0; seed < 5; seed++) {
        std::mt19937_64 random(1);
        std::uniform_real_distribution<float> dis(0, 1e9);
        std::uniform_real_distribution<double> dis2(0, 1e9);
        vector<float, AlignAlloc<float>> numbers(100000000);
        // vector<double, AlignAlloc<double>> numbers2(500000);
        for (float& e : numbers) e = dis(random);
        // for (double& e: numbers2) e = dis2(random);

        vector<function<double()>> func;
        /*func << L2(FSum<float>(numbers) / numbers.size());
        func << L2(NaiveSum<float>(numbers) / numbers.size());
        func << L2(ShardedSum2<float>(numbers) / numbers.size());
        func << L2(NaiveSum<float>(numbers) / numbers.size());
        func << L2(ShardedSum2<float>(numbers) / numbers.size());
        func << L2(ShardedSum<float>(numbers) / numbers.size());
        func << L2(NaiveSum<float>(numbers) / numbers.size());
        func << L2(NaiveSum<float>(numbers) / numbers.size());
        func << L2(ShardedSumSimdF(numbers) / numbers.size());
        func << L2(ShardedSumVectorizedF(numbers) / numbers.size());
        func << L2(NaiveSum<double>(numbers) / numbers.size());
        func << L2(NaiveSum<double>(numbers) / numbers.size());
        func << L2(ShardedSum<double>(numbers) / numbers.size());
        func << L2(ShardedSum2<double>(numbers) / numbers.size());
        func << L2(ShardedSum2<double>(numbers) / numbers.size());
        func << L2(NaiveSum<double>(numbers2) / numbers2.size());
        func << L2(ShardedSum2<double>(numbers2) / numbers2.size());
        func << L2(NaiveSum<double>(numbers2) / numbers2.size());
        func << L2(ShardedSum2<double>(numbers2) / numbers2.size());
        func << L2(NaiveSum<long double>(numbers) / numbers.size());
        func << L2(IncrementalMean<float>(numbers));
        func << L2(IncrementalMean<double>(numbers));*/

        func << L2(IncrementalMean<double>(numbers));
        func << L2(ShardedSum2<double>(numbers) / numbers.size());

        vector<string_view> name;
        /*name << "FSum<float>"sv;
        name << "NaiveSum<float>"sv;
        name << "ShardedSum2<float>"sv;
        name << "NaiveSum<float>"sv;
        name << "ShardedSum2<float>"sv;
        name << "ShardedSum<float>"sv;
        name << "NaiveSum<float>"sv;
        name << "NaiveSum<float>"sv;
        name << "ShardedSumSimdF"sv;
        name << "ShardedSumVectorizedF"sv;
        name << "NaiveSum<double>"sv;
        name << "NaiveSum<double>"sv;
        name << "ShardedSum<double>"sv;
        name << "ShardedSum2<double>"sv;
        name << "ShardedSum2<double>"sv;
        name << "xNaiveSum<double>"sv;
        name << "xShardedSum2<double>"sv;
        name << "xNaiveSum<double>"sv;
        name << "xShardedSum2<double>"sv;
        name << "NaiveSum<ldouble>"sv;
        name << "IncMean<float>"sv;
        name << "IncMean<double>"sv;*/

        name << "IncMean<double>"sv;
        name << "ShardedSum2<double>"sv;

        vector<ulong> ticks(name.size(), 0);
        vector<Accumulator<double>> acc(name.size());

        for (int cs = 0; cs < 10; cs++) {
            std::shuffle(numbers.begin(), numbers.end(), random);
            for (auto i : range(func.size())) {
                const auto& f = func[i];
                Timestamp begin;
                double result = f();
                Timestamp end;
                ticks[i] += begin.elapsed(end);
                acc[i] << result;
            }
        }

        ulong tmm = ticks[0];
        for (auto e : ticks) tmm = min<ulong>(tmm, e);
        double tm = tmm;

        vector<string> rows;
        for (auto i : range(func.size()))
            rows << format("%s|%.5f|%.5f|%f|%f", name[i], acc[i].mean(), acc[i].max() - acc[i].min(),
                           sqrt(acc[i].variance()), ticks[i] / tm);
        PrintTable(rows, '|', " ");
        println();
    }
}
