#pragma once
#include <core/std.h>
#include <sokoban/level.h>
#include <sokoban/state.h>

using Solution = vector<DynamicState>;

Solution Solve(const Level* level);

void GenerateDeadlocks(const Level* level);
