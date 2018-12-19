#pragma once
#include "int.h"
#include <functional>

constexpr double INF = std::numeric_limits<double>::infinity();

static_assert(sizeof(int) == 4);
static_assert(sizeof(long) == 8);
static_assert(sizeof(long long) == 8);
