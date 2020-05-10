#pragma once
#include <core/std.h>

enum class Figure : char { None, Dome, Player1, Player2 };

struct Cell {
    char level = 0;                // 2 bits
    Figure figure = Figure::None;  // 2 bits
};

inline bool operator==(Cell a, Cell b) { return a.level == b.level && a.figure == b.figure; }
inline bool operator!=(Cell a, Cell b) { return !(a == b); }

inline bool Less(Cell a, Cell b) {
    if (a.level != b.level) return a.level < b.level;
    return a.figure < b.figure;
}

inline size_t Hash(Cell a) { return int(a.figure) + (a.level << 2); }

inline Figure Other(Figure player) { return (player == Figure::Player1) ? Figure::Player2 : Figure::Player1; }

inline string_view PlayerName(Figure player) { return (player == Figure::Player1) ? "yellow" : "red"; }
