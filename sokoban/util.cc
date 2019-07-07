#include "sokoban/level.h"

inline bool free(Cell* a, const Boxes& boxes) {
	return a && !boxes[a->id];
}

// $$  $#  $#
// $$  $#  $$
static bool is_2x2_deadlock(Cell* box, const Boxes& boxes) {
	for (int d = 0; d < 4; d++) {
		Cell* a = box->dir[d];
		if (free(a, boxes))
			continue;
		Cell* b = box->dir[(d + 1) % 4];
		if (free(b, boxes))
			continue;
		if (!a && !b)
			return !box->goal;
		if (a) {
			Cell* c = a->dir[(d + 1) % 4];
			if (!free(c, boxes))
				return !(box->goal && a->goal && (!b || b->goal) && (!c || c->goal));
		}
		if (b) {
			Cell* c = b->dir[d];
			if (!free(c, boxes))
				return !(box->goal && b->goal && (!a || a->goal) && (!c || c->goal));
		}
	}
	return false;
}

// #$.
// .$#
static bool is_2x3_deadlock(Cell* pushed_box, const Boxes& boxes) {
    Cell* a = pushed_box;
    for (int d = 0; d < 4; d++) {
        Cell* b = a->dir[d];
        if (!b || !boxes[b->id])
            continue;
        if (a->goal && b->goal)
            continue;
        // Both A and B are boxes, and one of them is not on goal
        if (!a->dir[(d + 3) % 4] && !b->dir[(d + 1) % 4])
            return true;
        if (!a->dir[(d + 1) % 4] && !b->dir[(d + 3) % 4])
            return true;
    }
    return false;
}

bool is_simple_deadlock(Cell* pushed_box, const Boxes& boxes) {
	return is_2x2_deadlock(pushed_box, boxes) || is_2x3_deadlock(pushed_box, boxes);
}

bool is_frozen_on_goal(Cell* box, const Boxes& boxes) {
	if (!box->goal)
		return false;
	for (int d = 0; d < 4; d++) {
		Cell* a = box->dir[d];
		if (free(a, boxes))
			continue;
		Cell* b = box->dir[(d + 1) % 4];
		if (free(b, boxes))
			continue;
		if (!a && !b)
			return true;
		if (a) {
			Cell* c = a->dir[(d + 1) % 4];
			if (!free(c, boxes))
				return true;
		}
		if (b) {
			Cell* c = b->dir[d];
			if (!free(c, boxes))
				return true;
		}
	}
	return false;
}

static bool around(Cell* z, int side, const State& s, int s_dir) {
    Cell* m = z->dir[(s_dir + side) % 4];
    if (!m || s.boxes[m->id])
        return false;
    m = m->dir[s_dir];
    if (!m || s.boxes[m->id])
        return false;
    m = m->dir[s_dir];
    if (!m || s.boxes[m->id])
        return false;
    return true;
}

static bool around(Cell* z, const State& s, int s_dir) {
	return around(z, 1, s, s_dir) || around(z, 3, s, s_dir);
}

// can agent move to C without pushing any box?
bool is_cell_reachable(Cell* c, const State& s, small_bfs<Cell*>& visitor) {
	visitor.clear();
	visitor.add(c->level->cells[s.agent], s.agent);
    for (Cell* a : visitor)
		for (auto [_, b] : a->moves) {
			if (c == b)
				return true;
			if (!s.boxes[b->id])
				visitor.add(b, b->id);
		}
	return false;
}

bool is_reversible_push(const State& s, int dir, const Level* level, small_bfs<Cell*>& visitor) {
	Cell* agent = level->cells[s.agent];
	Cell* b = agent->dir[dir];
	Cell* c = b->dir[dir];
	if (!c || s.boxes[c->id])
		return false;

	if (around(agent, s, dir) || is_cell_reachable(c, s, visitor)) {
		State s2;
		s2.agent = b->id;
		s2.boxes = s.boxes;
		s2.boxes[b->id] = false;
		s2.boxes[s.agent] = true;

		Cell* b2 = b->dir[dir ^ 2];
		Cell* c2 = b2->dir[dir ^ 2];
		if (!c2 || s.boxes[c2->id])
			return false;

		return around(b, s2, dir ^ 2) || is_cell_reachable(c2, s2, visitor);
	}
	return false;
}

// frozen box deadlock
bool contains_frozen_boxes(Cell* agent, Boxes boxes, small_bfs<Cell*>& visitor) {
    int pushed_boxes = 0;
	int num_boxes = agent->level->num_goals;

    visitor.clear();
	visitor.add(agent, agent->id);

	for (Cell* a : visitor)
        for (auto [d, b] : a->moves) {
            if (visitor.visited[b->id])
                continue;
            if (!b->alive || !boxes[b->id]) {
                // agent moves to B
                visitor.add(b, b->id);
                continue;
            }

            Cell* c = b->dir[d];
            if (!c || !c->alive || boxes[c->id])
                continue;

			boxes[b->id] = false;
			boxes[c->id] = true;
            bool m = is_simple_deadlock(c, boxes);
            boxes[c->id] = false;
            if (m) {
				boxes[b->id] = true;
                continue;
            }

            // agent pushes box from B to C (and box disappears)
            if (--num_boxes == 1)
				return false;
            pushed_boxes += 1;
            visitor.clear();
			visitor.add(b, b->id);
            break;
        }

	// deadlock if any remaining box isn't on goal
	for (uint i = agent->level->num_goals; i < agent->level->num_alive; i++)
		if (boxes[i])
			return true;
	return false;
}

