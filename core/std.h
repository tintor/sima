#pragma once
#include <array>
#include <atomic>
#include <complex>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace fs = std::filesystem;

using std::complex;
using std::function;
using std::pair;
using std::regex;
using std::string;
using std::string_view;
using std::variant;

using std::array;

using std::map;
using std::set;
using std::unordered_map;
using std::unordered_set;

using std::multimap;
using std::multiset;
using std::unordered_multimap;
using std::unordered_multiset;

using std::cout;
using std::endl;
using std::istream;
using std::ostream;
using std::ostringstream;
using std::swap;

using std::numeric_limits;

using std::make_shared;
using std::nullopt;
using std::optional;
using std::shared_ptr;
using std::unique_ptr;

using std::atomic;

using namespace std::literals;

#include <core/span.h>
#include <core/vector.h>

#ifdef NDEBUG
constexpr bool debug = false;
#else
constexpr bool debug = true;
#endif

template <typename T>
void operator<<(vector<T>& v, const T& e) {
    v.push_back(e);
}

template <typename T, typename M>
void operator<<(vector<T>& v, const vector<M>& e) {
    v.insert(v.end(), e.begin(), e.end());
}

template <typename T, typename M>
void operator<<(vector<T>& v, span<M> e) {
    v.insert(v.end(), e.begin(), e.end());
}

// for std::visit and std::variant
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...)->overloaded<Ts...>;  // not needed as of C++20

constexpr double PI = M_PI;
constexpr double INF = numeric_limits<double>::infinity();

#include <core/callstack.h>

inline void Check(bool value, string_view message = "", const char* file = __builtin_FILE(),
                  unsigned line = __builtin_LINE(), const char* function = __builtin_FUNCTION()) {
    if (value) return;
    cout << "Check failed in " << function << " at " << file << ':' << line << " with message: " << message << endl;
    string s;
    Callstack().write(s, {});
    cout << s;
    cout.flush();
    exit(0);
}

#if 0
#define DCheck(A, B) Check(A, B)
#else
#define DCheck(A, B)
#endif

inline void Fail(string_view message = "", const char* file = __builtin_FILE(), unsigned line = __builtin_LINE(),
                 const char* function = __builtin_FUNCTION()) {
    cout << "Failed in " << function << " at " << file << ':' << line << " with message: " << message << endl;
    string s;
    Callstack().write(s, {});
    cout << s;
    cout.flush();
    exit(0);
}
