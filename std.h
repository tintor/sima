#pragma once
#include <string>
#include <string_view>
#include <utility>

#include <vector>
#include <array>

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

#include <iostream>

using std::string;
using std::string_view;
using std::pair;

using std::vector;
using std::array;

using std::unordered_map;
using std::unordered_set;
using std::map;
using std::set;

using std::cout;
using std::endl;

using namespace std::string_literals;

#ifdef NDEBUG
constexpr bool debug = false;
#else
constexpr bool debug = true;
#endif
