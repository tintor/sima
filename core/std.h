#pragma once
#include <string>
#include <string_view>
#include <regex>
#include <utility>
#include <complex>

#include <vector>
#include <array>

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

#include <iostream>
#include <limits>
#include <random>

using std::regex;
using std::string;
using std::string_view;
using std::pair;
using std::complex;

using std::vector;
using std::array;

using std::unordered_map;
using std::unordered_set;
using std::map;
using std::set;

using std::unordered_multimap;
using std::unordered_multiset;
using std::multimap;
using std::multiset;

using std::cout;
using std::endl;
using std::swap;

using std::numeric_limits;

using namespace std::literals;

#include <core/optional.h>
#include <core/span.h>

#ifdef NDEBUG
constexpr bool debug = false;
#else
constexpr bool debug = true;
#endif

constexpr double PI = M_PI;
constexpr double INF = numeric_limits<double>::infinity();

