#pragma once
#include <core/hash.h>
#include <core/std.h>
#include <santorini/cell.h>
#include <santorini/coord.h>
#include <absl/container/flat_hash_map.h>

#include <fstream>
#include <filesystem>

using Cells = array<Cell, 25>;

std::array<size_t, 7> zorbist[25];
struct InitZorbist {
    InitZorbist() {
        std::mt19937_64 re(0);
        FOR(i, 25) for (auto& e : zorbist[i]) e = re();
    }
} init_zorbist;

struct Board {
    bool setup = true;
    Figure player = Figure::Player1;
    optional<Coord> moved;
    bool built = false;

    Cells cell;

    const Cell& operator()(Coord c) const { return cell[c.v]; }
    Cell& operator()(Coord c) { return cell[c.v]; }

    size_t hash() const {
        if (hash_code != 0) return hash_code;

        size_t h = 0;
        for (int i = 0; i < 25; i++) {
            const auto& z = zorbist[i];
            auto c = cell[i];
            h ^= z[c.level];
            if (c.figure == Figure::Dome) h ^= z[4];
            if (c.figure == Figure::Player1) h ^= z[5];
            if (c.figure == Figure::Player2) h ^= z[6];
        }
        hash_code = h ? h : 1;
        return h;
    }

private:
    mutable size_t hash_code = 0;
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

void Render(const Board& board) {
    const std::string color[4] = { "\033[0m", "\033[0;34m", "\033[0;33m", "\033[1;31m" };
    const char* tower = ".#xo";

    for (int row = 0; row < 5; row++) {
        for (int col = 0; col < 5; col++) {
            const auto c = board.cell[row * 5 + col];
            // Color
            cout << color[c.level];
            // Figure
            if (c.figure == Figure::Dome) cout << '*';
            else if (c.figure == Figure::Player1) cout << 'A';
            else if (c.figure == Figure::Player2) cout << 'B';
            else if (c.figure == Figure::None) cout << tower[c.level];
            else cout << char(c.figure);
            cout << ' ';
        }
        cout << endl;
    }
    // Reset color
    cout << "\033[0m";
}

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

namespace std {
template <>
struct hash<Board> {
    size_t operator()(const Board& b) const { return b.hash(); }
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
    if (out.shape() != dim4(5, 5, 7)) Fail(string(out.shape()));
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

struct Score {
    uint p1 = 0, p2 = 0;
    Score() {}
    Score(const Score& s) : p1(s.p1), p2(s.p2) {}
    Score(uint p1, uint p2) : p1(p1), p2(p2) { }
    Score(Figure w) : p1((w == Figure::Player1) ? 1 : 0), p2((w == Figure::Player2) ? 1 : 0) {}

    float ValueP1() const { return p1 / float(p1 + p2); }
    Score operator+(Score o) { return {p1 + o.p1, p2 + o.p2}; }
    void operator+=(Score o) { *this = *this + o; }
};

class Values {
    static constexpr int Shards = 64;
   public:
    size_t Size() const { return Sum<size_t>(m_shard, [](const auto& e) { std::unique_lock lock(e._mutex); return e.data.size(); }); }

    /*const auto& Sample(std::mt19937_64& random) const {
        Check(m_data.size() > 0);
        auto it = (m_data.bucket_count() != m_last_buckets) ? m_data.begin() : m_last_iterator;
        std::uniform_int_distribution<int> dis(0, 9);
        for (auto i : range(dis(random))) {
            if (++it == m_data.end()) it = m_data.begin();
        }
        m_last_buckets = m_data.bucket_count();
        m_last_iterator = it;
        return *it;
    }*/

    void Clear() {
        for (auto& e : m_shard) {
            std::unique_lock lock(e._mutex);
            e.data.clear();
        }
    }

    void Add(const Board& board, Score score) {
        const Board b = Normalize(board);
        Shard& shard = m_shard[b.hash() % Shards];
        std::unique_lock lock(shard._mutex);
        shard.data.emplace(b, Score()).first->second += score;
    }

    void Merge(const Values& values) {
        FOR(i, Shards) {
            Shard& shard = m_shard[i];
            const Shard& o_shard = values.m_shard[i];
                
            std::unique_lock lock(shard._mutex);
            std::unique_lock lock2(o_shard._mutex);

            for (const auto& e : o_shard.data) {
                shard.data.emplace(e.first, Score()).first->second += e.second;
            }
        }
    }

    const Score* Lookup(const Board& board) const {
        const Board b = Normalize(board);
        const Shard& shard = m_shard[b.hash() % Shards];
        std::unique_lock lock(shard._mutex);
        auto it = shard.data.find(b);
        if (it == shard.data.end()) return nullptr;
        return &it->second;
    }

    float ValueP1(const Board& board) const {
        auto score = Lookup(board);
        return score ? score->ValueP1() : 0.5f;
    }

    void Export(std::filesystem::path filename) {
        vtensor out({5, 5, 7});
        std::ofstream os(filename);
        for (auto& shard : m_shard) {
            std::unique_lock lock(shard._mutex);
            for (const auto& [board, score] : shard.data) {
                ToTensor(board, out);
                FOR(x, 5) FOR(y, 5) {
                    tensor::type* s = out.data() + out.offset(x, y);
                    char b = 0;
                    FOR(i, 7) if (s[i] > 0) b |= 1 << i;
                    os.put(b);
                }
                os.write(reinterpret_cast<const char*>(&score.p1), sizeof(score.p1));
                os.write(reinterpret_cast<const char*>(&score.p2), sizeof(score.p2));
            }
        }
    }

   private:
    using Data = absl::flat_hash_map<Board, Score>;
    struct Shard {
        Data data;
        mutable mutex _mutex;
    } m_shard[Shards];
    
    //mutable size_t m_last_buckets = 0;
    //mutable Data::const_iterator m_last_iterator;

};
