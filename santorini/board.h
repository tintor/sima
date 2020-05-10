#pragma once
#include <core/hash.h>
#include <core/std.h>
#include <santorini/cell.h>
#include <santorini/coord.h>

struct Board {
    bool setup = true;
    Figure player = Figure::Player1;
    optional<Coord> moved;
    bool built = false;

    array<Cell, 25> cell;

    const Cell& operator()(Coord c) const { return cell[c.v]; }
    Cell& operator()(Coord c) { return cell[c.v]; }
};

// Symmetrical boards are equal!
bool Equal(const Board& a, const Board& b) {
    if (a.setup != b.setup || a.player != b.player || a.moved != b.moved || a.built != b.built) return false;
    for (int t = 0; t < 8; t++) {
        if (All(kAll, [&](Coord e) { return a.cell[e.v] == b.cell[Transform(e, t).v]; })) return true;
    }
    return false;
}

size_t Hash(const array<Cell, 25>& cell, int transform) {
    hash h;
    for (Coord e : kAll) h << Hash(cell[Transform(e, transform).v]);
    return h.seed;
}

// Symmetrical boards will have the same hash!
size_t Hash(const Board& a) {
    hash h;
    h << a.setup << int(a.player) << a.moved.has_value();
    if (a.moved.has_value()) h << a.moved->v;
    h << a.built;

    size_t m = 0;
    for (int t = 0; t < 8; t++) m ^= Hash(a.cell, t);
    h << m;
    return h.seed;
}

template <typename Fn>
int Count(const Board& board, const Fn& fn) {
    return CountIf(kAll, L(fn(board(e))));
}

// Network input description:
// Board:
// 0/1 - setup
// 0/1 - moved
// 0/1 - built
// 0/1 - player
// 25x - cell
// Cell:
// 0/3 - level
// 0/1 - dome
// 0/1 - player1
// 0/1 - player2

// Output:
// 0..1 - value function
