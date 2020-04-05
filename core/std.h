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
#include <optional>

#include <atomic>

using std::regex;
using std::string;
using std::string_view;
using std::pair;
using std::complex;

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

using std::optional;
using std::nullopt;
using std::unique_ptr;

using std::atomic;

using namespace std::literals;

// Define vector with sizeof(vector) == 16 instead of 24
#include <boost/container/vector.hpp>
template<typename T, typename Allocator = std::allocator<T>>
using vector = boost::container::vector<T, Allocator, typename boost::container::vector_options<boost::container::stored_size<uint32_t>>::type>;
static_assert(sizeof(vector<int>) == 16);

// vector implemented in array
#include <boost/container/static_vector.hpp>
using boost::container::static_vector;

// static storage for small number of elements
#include <boost/container/small_vector.hpp>
using boost::container::small_vector;

#include <core/span.h>

#ifdef NDEBUG
constexpr bool debug = false;
#else
constexpr bool debug = true;
#endif

template<typename T>
void operator<<(vector<T>& v, const T& e) {
	v.push_back(e);
}

template<typename T, typename M>
void operator<<(vector<T>& v, const vector<M>& e) {
	v.insert(v.end(), e.begin(), e.end());
}

template<typename T, typename M>
void operator<<(vector<T>& v, span<M> e) {
	v.insert(v.end(), e.begin(), e.end());
}

constexpr double PI = M_PI;
constexpr double INF = numeric_limits<double>::infinity();

