#pragma once
#include <core/std.h>

struct Coord {
    char v;

    constexpr Coord() : v(-1) {}
    constexpr Coord(int x, int y) : v((x >= 0 && y >= 0 && x < 5 && y < 5) ? y * 5 + x : -1) {}

    Coord(const Coord& e) : v(e.v) {}
    Coord(Coord&& e) : v(e.v) {}

    void operator=(const Coord& e) { v = e.v; }
    void operator=(Coord&& e) { v = e.v; }

    int x() const { return v % 5; }
    int y() const { return v / 5; }
};

inline bool operator==(Coord a, Coord b) { return a.v == b.v; }
inline bool operator!=(Coord a, Coord b) { return !(a == b); }

inline Coord Transpose(Coord e) { return {e.y(), e.x()}; }

inline Coord FlipX(Coord e) { return {4 - e.x(), e.y()}; }

inline Coord FlipY(Coord e) { return {e.x(), 4 - e.y()}; }

inline Coord Transform(Coord e, int code) {
    if (code & 1) e = FlipX(e);
    if (code & 2) e = FlipY(e);
    if (code & 4) e = Transpose(e);
    return e;
}

inline bool IsValid(Coord a) { return a.v != -1; }

inline constexpr bool Nearby(Coord src, Coord dest) {
    int x = src.v % 5 - dest.v % 5;
    int y = src.v / 5 - dest.v / 5;
    return -1 <= x && x <= 1 && -1 <= y && y <= 1;
}

static_assert(Nearby({1, 1}, {2, 2}));
static_assert(Nearby({2, 1}, {1, 2}));

inline vector<Coord> All() {
    vector<Coord> out;
    out.reserve(25);
    for (int x = 0; x < 5; x++)
        for (int y = 0; y < 5; y++) out.emplace_back(x, y);
    return out;
}

const vector<Coord> kAll = All();
