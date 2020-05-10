#pragma once
#include <array>
#include <atomic>
#include <complex>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <regex>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using std::complex;
using std::pair;
using std::regex;
using std::string;
using std::string_view;

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
using std::swap;

using std::numeric_limits;

using std::nullopt;
using std::optional;
using std::unique_ptr;

using std::atomic;

using namespace std::literals;

#include <core/vector.h>
#include <core/span.h>

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

constexpr double PI = M_PI;
constexpr double INF = numeric_limits<double>::infinity();
