#pragma once
#include <core/hash.h>
#include <core/std.h>
#include <santorini/cell.h>
#include <santorini/coord.h>
#include <absl/container/flat_hash_map.h>

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

ostream& operator<<(ostream& os, const Board& board) {
    for (int row = 0; row < 5; row++) {
        for (int col = 0; col < 5; col++) {
            const Cell c = board.cell[row * 5 + col];
            os << char(c.figure) << char(c.level ? '0' + c.level : '.') << ' ';
        }
        os << endl;
    }
    char p[2] = {char(board.player), 0};
    os << format("setup %s, player %s", board.setup, p);
    if (board.moved) os << format(" moved %s%s", board.moved->x(), board.moved->y());
    os << format(" built %s", board.built) << endl;
    return os;
}

void Print(const Board& board) { cout << board; }


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

void Serialize(ostream& os, Cell cell) {
    os << Bit(cell.level == 0);
    os << Bit(cell.level == 1);
    os << Bit(cell.level == 2);
    os << Bit(cell.level == 3);
    os << Bit(cell.figure == Figure::Dome);
    os << Bit(cell.figure == Figure::Player1);
    os << Bit(cell.figure == Figure::Player2);
}

void Serialize(ostream& os, const Board& board) {
    os << Bit(board.setup);
    for (Coord e : kAll) os << Bit(board.moved && *board.moved == e);
    os << Bit(board.built);
    os << Bit(board.player == Figure::Player1);
    for (Coord e : kAll) Serialize(os, board(e));
}

inline void Serialize(const Board& board, tensor out) {
    Check(out.rank == 1);
    Check(out.size == BoardBits);
    ostringstream os;
    os.str().reserve(BoardBits);
    Serialize(os, board);
    Check(os.str().size() == BoardBits);
    for (size_t i = 0; i < BoardBits; i++) out(i) = (os.str()[i] == '1') ? 1.f : 0.f;
}

class Values {
   public:
    struct Wins {
        uint32_t player1, player2;
        float Value() const { return player1 / double(player1 + player2); }
    };

    size_t Size() const { return m_data.size(); }

    const auto& Sample(std::mt19937_64& random) const {
        Check(m_data.size() > 0);
        auto it = (m_data.bucket_count() != m_last_buckets) ? m_data.begin() : m_last_iterator;
        std::uniform_int_distribution<int> dis(0, 9);
        for (auto i : range(dis(random))) {
            if (++it == m_data.end()) it = m_data.begin();
        }
        m_last_buckets = m_data.bucket_count();
        m_last_iterator = it;
        return *it;
    }

    void Add(const Board& board, uint32_t wins1, uint32_t wins2) {
        auto result = m_data.emplace(Normalize(board), Wins{0, 0});
        auto& wins = result.first->second;
        wins.player1 += wins1;
        wins.player2 += wins2;
    }

    const Wins* Lookup(const Board& board) const {
        auto it = m_data.find(Normalize(board));
        if (it == m_data.end()) return nullptr;
        return &it->second;
    }

    double Value(const Board& board) const {
        auto it = m_data.find(Normalize(board));
        if (it == m_data.end()) return 0.5;
        return it->second.Value();
    }

    void Export(string_view filename, int wmin) {
        // std::ofstream os((string(filename)));
        vector<Board> outs;
        for (const auto& [board, wins] : m_data) {
            if (wins.player1 + wins.player2 < wmin) continue;
            outs.clear();
            double w = wins.player1 / double(wins.player1 + wins.player2);
            for (int transform = 0; transform < 8; transform++) {
                // TODO(Marko) deduplicate
                Board out;
                out.cell = Transform(board.cell, transform);
                if (!contains(outs, out)) outs << out;
            }
            for (const auto& out : outs) {
                Print(out);
                println("wins: %s %s", wins.player1, wins.player2);
                // os << out << ' ' << std::setprecision(10) << w << '\n';
            }
        }
    }

   private:
    using Data = absl::flat_hash_map<Board, Wins>;
    Data m_data;
    mutable size_t m_last_buckets = 0;
    mutable Data::const_iterator m_last_iterator;

};
