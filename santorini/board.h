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
// Board (5x5xCell):
// Cell:
// 4 - level
// 1 - dome
// 1 - ego
// 1 - opponent

constexpr int CellBits = 4 + 3;
constexpr int BoardBits = 25 * CellBits;

void ToTensor(const Board& board, tensor out) {
    Check(out.shape() == dim4(5, 5, 7));
    const auto other = Other(board.player);
    FOR(x, 5) FOR(y, 5) {
        tensor::type* s = out.data() + out.offset(x, y);
        auto cell = board.cell[y * 5 + x];
        s[0] = cell.level == 0;
        s[1] = cell.level == 1;
        s[2] = cell.level == 2;
        s[3] = cell.level == 3;
        s[4] = cell.figure == Figure::Dome;
        s[5] = cell.figure == other;
        s[6] = cell.figure == board.player;
    }
}

class Values {
   public:
    struct Score {
        uint32_t wins, games;
        float Value() const { return wins / float(games); }
    };

    size_t Size() const { return m_data.size(); }

    auto begin() const { return m_data.begin(); }
    auto end() const { return m_data.end(); }

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

    void Add(const Board& board, uint32_t wins, uint32_t games) {
        auto result = m_data.emplace(Normalize(board), Score{0, 0});
        auto& score = result.first->second;
        score.wins += wins;
        score.games += games;
    }

    const Score* Lookup(const Board& board) const {
        auto it = m_data.find(Normalize(board));
        if (it == m_data.end()) return nullptr;
        return &it->second;
    }

    float Value(const Board& board) const {
        auto it = m_data.find(Normalize(board));
        if (it == m_data.end()) return 0.5f;
        return it->second.Value();
    }

    void Export(string_view filename) {
        vtensor out({5, 5, 7});
        std::ofstream os(filename);
        for (const auto& [board, score] : m_data) {
            ToTensor(board, out);
            FOR(x, 5) FOR(y, 5) {
                tensor::type* s = out.data() + out.offset(x, y);
                char b = 0;
                FOR(i, 7) if (s[i] > 0) b |= 1 << i;
                os.put(b);
            }
            os.write(reinterpret_cast<const char*>(&score.wins), sizeof(score.wins));
            os.write(reinterpret_cast<const char*>(&score.games), sizeof(score.games));
        }
    }

   private:
    using Data = absl::flat_hash_map<Board, Score>;
    Data m_data;
    mutable size_t m_last_buckets = 0;
    mutable Data::const_iterator m_last_iterator;

};
