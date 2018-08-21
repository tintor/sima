#pragma once
#include <string>
#include <string_view>
#include <regex>
#include <utility>

#include <vector>
#include <array>

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

#include <iostream>

using std::regex;
using std::string;
using std::string_view;
using std::pair;

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

using namespace std::literals;

#ifdef NDEBUG
constexpr bool debug = false;
#else
constexpr bool debug = true;
#endif

template<typename T>
class optional {
public:
	optional() : _has_value(false) { }
	optional(T value) : _has_value(true), _value(value) { }
	bool has_value() const { return _has_value; }
	T& operator*() { assert(_has_value); return _value; }
	T* operator->() { assert(_has_value); return &_value; }
private:
	bool _has_value;
	T _value;
};
