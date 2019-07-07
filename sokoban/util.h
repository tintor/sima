#include "sokoban/level.h"

bool is_simple_deadlock(Cell* box, const Boxes& boxes);
bool is_frozen_on_goal(Cell* box, const Boxes& boxes);
bool is_reversible_push(const State& s, int s_dir, const Level* level, small_bfs<Cell*>& visitor);
bool is_cell_reachable(Cell* c, const State& s, small_bfs<Cell*>& visitor);
bool contains_frozen_boxes(Cell* agent, Boxes boxes, small_bfs<Cell*>& visitor);
