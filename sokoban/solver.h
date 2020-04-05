#pragma once
#include <core/std.h>
#include <sokoban/state.h>
#include <sokoban/level.h>

using Solution = vector<DynamicState>;

Solution Solve(const Level* level);
