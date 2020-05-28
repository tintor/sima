#pragma once
#include <core/hash.h>
#include <core/std.h>
#include <core/arrayfire.h>
#include <santorini/cell.h>
#include <santorini/coord.h>

#include <fstream>

using Cells = array<Cell, 25>;

struct Board {
    bool setup = true;
    Figure player = Figure::Player1;
    optional<Coord> moved;
    bool built = false;

    Cells cell;

    const Cell& operator()(Coord c) const { return cell[c.v]; }
    Cell& operator()(Coord c) { return cell[c.v]; }
};

bool Less(const Cells& a, const Cells& b) {
    for (int i = 0; i < 25; i++) {
        if (a[i] != b[i]) return Less(a[i], b[i]);
    }
    return false;
}

Cells Transform(const Cells& cell, int transform) {
    Cells out;
    for (Coord e : kAll) out[e.v] = cell[Transform(e, transform).v];
    return out;
}

// Returns the same board for all symmetry variations!
Board Normalize(const Board& board) {
    // generate all 8 transforms and return the smallest one
    Board out = board;
    for (int transform = 1; transform < 8; transform++) {
        Cells m = Transform(board.cell, transform);
        if (Less(m, out.cell)) out.cell = m;
    }
    return out;
}

bool operator==(const Board& a, const Board& b) {
    return a.setup == b.setup && a.player == b.player && a.moved == b.moved && a.built == b.built && a.cell == b.cell;
}

size_t Hash(const Cells& cell) {
    hash h;
    for (Coord e : kAll) h << Hash(cell[e.v]);
    return h.seed;
}

size_t Hash(const Board& a) {
    hash h;
    h << a.setup << int(a.player) << a.moved.has_value();
    if (a.moved.has_value()) h << a.moved->v;
    h << a.built << Hash(a.cell);
    return h.seed;
}

namespace std {
template <>
struct hash<Board> {
    size_t operator()(const Board& b) const { return Hash(b); }
};
}  // namespace std

template <typename Fn>
int Count(const Board& board, const Fn& fn) {
    return CountIf(kAll, L(fn(board(e))));
}

// Network input description:
// Board (28 + 25xCell):
// 1 - setup
// 25 - which figure moved (if any)
// 1 - built
// 1 - player
// 25x - cell

// Cell (7):
// 4 - level
// 1 - dome
// 1 - player1
// 1 - player2

constexpr int CellBits = 4 + 3;
constexpr int BoardBits = 3 + 25 * (1 + CellBits);

using std::ostream;

inline char Bit(bool a) { return a ? '1' : '0'; }

ostream& operator<<(ostream& os, Cell cell) {
    os << Bit(cell.level == 0);
    os << Bit(cell.level == 1);
    os << Bit(cell.level == 2);
    os << Bit(cell.level == 3);
    os << Bit(cell.figure == Figure::Dome);
    os << Bit(cell.figure == Figure::Player1);
    os << Bit(cell.figure == Figure::Player2);
    return os;
}

ostream& operator<<(ostream& os, const Board& board) {
    os << Bit(board.setup);
    for (Coord e : kAll) os << Bit(board.moved && *board.moved == e);
    os << Bit(board.built);
    os << Bit(board.player == Figure::Player1);
    for (Coord e : kAll) os << board(e);
    return os;
}

void Serialize(const Board& board, af::array& out) {
    Check(out.dims() == 1);
    Check(out.elements() == BoardBits);
    ostringstream os;
    os.str().reserve(BoardBits);
    os << board;
    Check(os.str().size() == BoardBits);
    for (size_t i = 0; i < BoardBits; i++) out(i) = (os.str()[i] == '1') ? 1.f : 0.f;
}

class Values {
   private:
    struct Wins {
        uint32_t player1, player2;
    };

   public:
    size_t Size() const { return m_data.size(); }

    void Add(const Board& board, uint32_t wins1, uint32_t wins2) {
        auto result = m_data.emplace(Normalize(board), Wins{0, 0});
        auto& wins = result.first->second;
        wins.player1 += wins1;
        wins.player2 += wins2;
    }

    optional<double> Value(const Board& board) const {
        auto it = m_data.find(Normalize(board));
        if (it == m_data.end()) return nullopt;
        const auto& wins = it->second;
        return ((board.player == Figure::Player1) ? wins.player1 : wins.player2) / double(wins.player1 + wins.player2);
    }

    void Export(string_view filename) {
        std::ofstream os((string(filename)));
        for (const auto& [board, wins] : m_data) {
            double w = wins.player1 / double(wins.player1 + wins.player2);
            for (int transform = 0; transform < 8; transform++) {
                // TODO(Marko) deduplicate
                Board out;
                out.cell = Transform(board.cell, transform);
                os << out << ' ' << std::setprecision(10) << w << '\n';
            }
        }
    }

   private:
    unordered_map<Board, Wins> m_data;
};
