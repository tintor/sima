#pragma once
#include <core/std.h>

#include <algorithm>

inline double env(const char* name, double def) {
    const char* v = std::getenv(name);
    return v ? std::atof(v) : def;
}

template<typename T>
T RoundUp(T a, T b) {
    auto m = a % b;
    return m ? a + b - m : a;
}

#define LAMBDA(EXPR) \
    (const auto& e) { return EXPR; }

#define L(A) [&](const auto& e) { return A; }

template <typename T>
int sign(T a) {
    return (0 < a) - (a < 0);
}

template <typename T>
T sqr(T a) {
    return a * a;
}

template <typename T>
T cube(T a) {
    return a * a * a;
}

template <typename T, typename P>
T lerp(T a, T b, P t) {
    return a + (b - a) * t;
}

template <size_t S>
bool aligned(const void* ptr) {
    return (reinterpret_cast<size_t>(ptr) % S) == 0;
}

template <typename T>
inline bool ordered(T a, T b, T c) {
    return a <= b && b <= c;
}

template <typename T>
void minimize(T& a, T b) {
    if (b < a) a = b;
}

template <typename T>
void maximize(T& a, T b) {
    if (b > a) a = b;
}

template <typename T>
T min(T a, T b) {
    return (a < b) ? a : b;
}

template <typename T, typename... Args>
T min(T a, T b, Args... args) {
    return min(a, min(b, args...));
}

template <typename T>
T max(T a, T b) {
    return (a > b) ? a : b;
}

template <typename T, typename... Args>
T max(T a, T b, Args... args) {
    return max(a, max(b, args...));
}

template <typename Container>
auto min(const Container& container) {
    return *std::min_element(container.begin(), container.end());
}

template <typename Container>
auto max(const Container& container) {
    return *std::max_element(container.begin(), container.end());
}

template <typename T>
T clamp(T a, T low, T high) {
    if (a < low) return low;
    if (a > high) return high;
    return a;
}

template <typename T>
T median(T x, T y, T z) {
    if (x <= y) {
        if (z <= x) return x;
        if (y <= z) return y;
        return z;
    }
    // y < x
    if (z <= y) return y;
    if (x <= z) return x;
    return z;
}

using std::sort;

template <typename Container>
void sort(Container& container) {
    std::sort(container.begin(), container.end());
}

template <typename Container, typename Func>
void sort(Container& container, const Func& func) {
    std::sort(container.begin(), container.end(), func);
}

template <typename Vector>
void remove_dups(Vector& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

template <typename Vector, typename Func>
void remove_dups(Vector& v, const Func& less, const Func& equal) {
    std::sort(v.begin(), v.end(), less);
    v.erase(std::unique(v.begin(), v.end(), equal), v.end());
}

template <typename T>
bool contains(const vector<T>& container, T element) {
    for (const T& e : container)
        if (e == element) return true;
    return false;
}

template <typename T>
bool contains(cspan<T> container, T element) {
    for (const T& e : container)
        if (e == element) return true;
    return false;
}

template <typename T, size_t M>
bool contains(const array<T, M>& container, T element) {
    return contains(cspan<T>(container), element);
}

template <typename Map, typename Value>
inline bool contains_value(const Map& map, const Value& value) {
    for (const auto& [k, v] : map)
        if (v == value) return true;
    return false;
}

template <typename Iterable, typename UnaryPredicate>
bool All(const Iterable& iterable, const UnaryPredicate& predicate) {
    for (const auto& e : iterable)
        if (!predicate(e)) return false;
    return true;
}

template <typename Iterable, typename UnaryPredicate>
bool Any(const Iterable& iterable, const UnaryPredicate& predicate) {
    for (const auto& e : iterable)
        if (predicate(e)) return true;
    return false;
}

template <typename Iterable, typename UnaryPredicate>
bool None(const Iterable& iterable, const UnaryPredicate& predicate) {
    for (const auto& e : iterable)
        if (predicate(e)) return false;
    return true;
}

template <typename Iterable, typename Value, typename Accumulate>
Value Accum(const Iterable& iterable, const Value& init, const Accumulate& accumulate) {
    Value acc = init;
    for (const auto& e : iterable) accumulate(acc, e);
    return acc;
}

// Welford algorithm for computing variance.
// More numerically stable than naive sum of squares variance.
template <typename T>
class Accumulator {
   public:
    void operator<<(T a) {
        m_count += 1;
        double delta = a - m_mean;
        m_mean += delta / m_count;
        double delta2 = a - m_mean;
        m_m2 += delta * delta2;

        if (a < m_min) m_min = a;
        if (a > m_max) m_max = a;
    }

    T mean() const { return m_mean; }
    T variance() const { return m_m2 / m_count; }
    T stdev() const { return sqrt(m_m2 / m_count); }
    T min() const { return m_min; }
    T max() const { return m_max; }

   private:
    size_t m_count = 0;
    double m_mean = 0;
    double m_m2 = 0;
    T m_min = std::numeric_limits<T>::infinity();
    T m_max = -std::numeric_limits<T>::infinity();
};

template <typename T>
struct Aggregates {
    Aggregates(const T* begin, const T* end) {
        double sum = 0;
        T rmin = std::numeric_limits<T>::infinity();
        T rmax = -std::numeric_limits<T>::infinity();
        for (auto p = begin; p != end; p++) {
            sum += *p;
            if (*p < rmin) rmin = *p;
            if (*p > rmax) rmax = *p;
        }
        mean = sum / (end - begin);
        min = rmin;
        max = rmax;

        double m2 = 0;
        for (auto p = begin; p != end; p++) m2 += sqr(*p - mean);
        variance = m2 / (end - begin);
    }

    double mean;
    double variance;
    T min;
    T max;
};

template <typename Iterable, typename UnaryPredicate>
size_t CountIf(const Iterable& iterable, const UnaryPredicate& predicate) {
    size_t count = 0;
    for (const auto& e : iterable)
        if (predicate(e)) count += 1;
    return count;
}

template <typename T, typename P>
size_t IndexOfMax(cspan<T> s, P measure) {
    auto mv = measure(s[0]);
    int mi = 0;
    for (size_t i = 1; i < s.size(); i++) {
        auto v = measure(s[i]);
        if (v > mv) {
            mv = v;
            mi = i;
        }
    }
    return mi;
}

template <typename T, typename Range>
T Product(const Range& range) {
    T out = 1;
    for (const auto& e : range) out *= e;
    return out;
}

template <typename T, typename Range>
T Sum(const Range& range) {
    T out = 0;
    for (const auto& e : range) out += e;
    return out;
}

template <typename CollectionA, typename CollectionB>
void Copy(const CollectionA& a, CollectionB& b) {
    Check(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) b[i] = a[i];
}

string Demangle(const char* name);

template <typename T>
inline string TypeName(const T& value) {
    return Demangle(typeid(value).name());
}

template <typename Type, typename Base>
inline bool IsType(shared_ptr<Base> a) {
    const auto& e = *a.get();
    return typeid(e) == typeid(Type);
}

void PrintTable(cspan<string> rows, char separator, string_view splitter, cspan<int> right = {});
