#pragma once
#include <core/int.h>
#include <core/std.h>

struct Level;

struct Cell {
    const Level* level;
    int id;
    int xy;
    bool goal;
    bool sink;
    bool alive;

    array<Cell*, 4> _dir;
    Cell* dir(int d) const { return _dir[d & 3]; }

    vector<pair<int, Cell*>> moves;
    vector<pair<Cell*, Cell*>> pushes;  // (box_dest, agent_src)

    array<Cell*, 8> dir8;

    constexpr static uint Inf = numeric_limits<int>::max();
    vector<uint> push_distance;  // from any alive cell to any goal cell (or Inf if not reachable)
    uint min_push_distance;      // min(push_distance)

    bool straight() const { return moves.size() == 2 && (moves[0].first ^ 2) == moves[1].first; }
};
