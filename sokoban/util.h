#pragma once
#include "sokoban/level.h"

bool is_simple_deadlock(Cell* box, const Boxes& boxes);

// Assumes no deadlocks!
bool is_frozen_on_goal_simple(Cell* box, const Boxes& boxes);
Boxes goals_with_frozen_boxes(Cell* agent, const Boxes& boxes, small_bfs<Cell*>* visitor = nullptr);

bool is_reversible_push(const State& s, int s_dir, const Level* level, small_bfs<Cell*>& visitor);
bool is_cell_reachable(Cell* c, const State& s, small_bfs<Cell*>& visitor);
bool contains_frozen_boxes(Cell* agent, Boxes boxes, small_bfs<Cell*>& visitor);
