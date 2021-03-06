#define CATCH_CONFIG_RUNNER
#include <core/dynamic_array.h>
#include <core/each.h>
#include <core/exception.h>
#include <core/format.h>
#include <core/small_bfs.h>
#include <core/util.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <catch.hpp>
#include <magic_enum.hpp>

enum class God : char {
    Dead,
    None,
    // Simple
    Apollo,
    Artemis,
    Athena,
    Atlas,
    Demeter,
    Hephaestus,
    /*partial*/ Hermes,
    Minotaur,
    Pan,
    Prometheus,
    // Advanced
    Ares,
    Bia,
    Selene,
    Eros,
    Chronus,
    Hera,
    Limus,
    Medusa,
    /*partial*/ Morpheus,
    Poseidon,
    /*partial*/ Triton,
    Zeus,
    // - Aphrodite!, Chaos!, Charon*, Circe!, Dionysus*, Hestia*, Hypnus*, Persephone,
    // Golden Fleece
    Hades,
    // - Aeolus*, Charybdis*, Clio*, EuropaTalus*, Gaea!, Graeae*, Harpies*, Hecate, Moerae*, Nemesis, Siren*,
    // - Tartarus*, Terpsichore, Urania*
    // Heroes
    // - Achilles*, Adonis*, Atalanta, Bellerophon*, Heracles*, Jason*, Medea*, Odysseus, Polyphemus*, Theseus*
    // Underworld
    Pegasus
};
// - Tyche, Sculla*, Proteus*, CastorPollux*, Eris, Maenads*, Asteria*, Hippolyta*, Hydra*, Iris*, Nyx

using Coord = char;  // 0-24 or -1 if empty

struct Cell {
    char level = 0;  // 0-3 height of level

    char figure = ' ';
    // <space> - empty
    // ) - dome
    // A - first player, B - second player, ... (uppercase male, lowercase female)
    // @ - golden fleece
    // T - talus

    // AddToken only:
    // $ - coin
    // * - whirlpool
    // T - talus
    // ~ - wind
    // > - siren's arrow
    // O - abyss
    // # - fate

    bool whirlpool = false;
    bool coin = false;
    bool abyss = false;
    bool fate = false;
};

bool Dome(Cell c) { return c.figure == ')'; }

bool Empty(Cell c) { return c.figure == ' '; }

bool operator<(Cell a, Cell b) {
    if (a.level != b.level) return a.level < b.level;
    return a.figure < b.figure;
}

bool operator==(Cell a, Cell b) { return a.level == b.level && a.figure == b.figure; }

bool operator!=(Cell a, Cell b) { return !(a == b); }

struct State {
    array<array<Cell, 5>, 5> cell;
    array<God, 2> gods;
    char player = 0;  // player about to play
    bool victory = false;
    bool setup = true;
    bool athenaMovedUp = false;
    God golden_fleece = God::None;
    char morpheus_blocks = 0;
    char wind = -1;
    char arrow = -1;

    // transient
    char ares_removed_blocks = 0;
    bool hermes_slow_move = false;
    Coord lastMoveSrc = -1;
    Coord lastMove = -1;   // coordinate of builder who moved last
    Coord lastBuild = -1;  // coordinate of last build
    char builds = 0;
    char moves = 0;

    Cell operator[](int a) const { return cell[a / 5][a % 5]; }
    Cell& operator[](int a) { return cell[a / 5][a % 5]; }

    State() {}
    State(string_view str) {
        for (int i = 0; i < 25; i++) {
            char l = str[i * 2];
            REQUIRE((l == '.' || ('0' <= l && l <= '3')));
            (*this)[i].level = (l == '.') ? 0 : (l - '0');

            char f = str[i * 2 + 1];
            REQUIRE((f == ' ' || f == ')' || ('a' <= f && f <= 'd') || ('A' <= f && f <= 'D')));
            (*this)[i].figure = f;
        }
    }
};

auto Cells() {
    array<int, 25> data;
    for (int i = 0; i < data.size(); i++) data[i] = i;
    return data;
}

bool IsMale(char figure) { return 'A' <= figure && figure <= 'D'; }
bool IsMale(Cell c) { return IsMale(c.figure); }

bool IsFemale(char figure) { return 'a' <= figure && figure <= 'd'; }
bool IsFemale(Cell c) { return IsFemale(c.figure); }

int Player(char figure) {
    if (IsMale(figure)) return figure - 'A';
    if (IsFemale(figure)) return figure - 'a';
    return -1;
}

int Player(Cell c) { return c.figure; }

int PrevPlayer(const State& state) {
    for (int i = 1; i < state.gods.size(); i++) {
        int p = (state.player + state.gods.size() - i) % state.gods.size();
        if (state.gods[p] != God::Dead) return p;
    }
    THROW(runtime_error);
}

int NextPlayer(const State& state) {
    for (int i = 1; i < state.gods.size(); i++) {
        int p = (state.player + i) % state.gods.size();
        if (state.gods[p] != God::Dead) return p;
    }
    THROW(runtime_error);
}

auto PrevPlayerCells(const State& state) {
    int p = PrevPlayer(state);
    static_vector<int, 4> data;
    for (int i = 0; i < 25; i++)
        if (Player(state[i]) == p) data.push_back(i);
    return data;
}

auto PlayerCells(const State& state) {
    static_vector<int, 4> data;
    for (int i = 0; i < 25; i++)
        if (Player(state[i]) == state.player) data.push_back(i);
    return data;
}

// TODO precompute
auto CellsAround(Coord a, bool includeCenter = false) {
    static_vector<Coord, 25> data;
    int row_min = max(0, a / 5 - 1), row_max = min(4, a / 5 + 1);
    int col_min = max(0, a % 5 - 1), col_max = min(4, a % 5 + 1);
    for (int r = row_min; r <= row_max; r++)
        for (int c = col_min; c <= col_max; c++) {
            Coord b = r * 5 + c;
            if (b != a || includeCenter) data.push_back(b);
        }
    return data;
}

// TODO precompute
auto CellsAroundL2(Coord a) {
    static_vector<Coord, 25> data;
    int row = a / 5;
    int col = a % 5;

    array<pair<int, int>, 8> jumps = {
        pair<int, int>{1, 2}, {1, -2}, {-1, 2}, {-1, -2}, {2, 1}, {-2, 1}, {2, -1}, {-2, -1},
    };
    for (auto p : jumps) {
        int r = row + p.first;
        int c = col + p.second;
        if (r >= 0 && r <= 4 && c >= 0 && c <= 4) data.push_back(r * 5 + c);
    }
    return data;
}

// TODO precompute
Coord CellInDirection(Coord a, Coord b) {
    int row = b / 5 + b / 5 - a / 5;
    int col = b % 5 + b % 5 - a % 5;
    return (0 <= row && row < 5 && 0 <= col && col < 5) ? row * 5 + col : -1;
}

// TODO precompute
bool OnPerimeter(Coord a) {
    int row = a / 5;
    int col = a % 5;
    return row == 0 || row == 4 || col == 0 || col == 4;
}

// NON Athena, NON Minotaur, NON Pan
State MoveBuilder(const State& state, Coord ia, Coord ib) {
    State s = state;
    s.lastMove = ib;
    if (state[ib].level == 3 && state[ia].level < 3) {
        s.victory = true;
        // Hera prevents others from winning by climbing on perimeter tower
        if (state.gods[state.player] != God::Hera && contains(state.gods, God::Hera) && !OnPerimeter(ib))
            s.victory = false;
    }
    swap(s[ia].figure, s[ib].figure);
    return s;
}

Coord LMiddle1(Coord a, int dr, int dc) {
    if (dr == 2) return a + (dr - 1) * 5 + dc;
    if (dr == -2) return a + (dr + 1) * 5 + dc;
    if (dc == 2) return a + dr * 5 + dc - 1;
    if (dc == -2) return a + dr * 5 + dc + 1;
    return -1;
}

Coord LMiddle2(Coord a, int dr, int dc) {
    if (dr > 0)
        dr -= 1;
    else
        dr += 1;
    if (dc > 0)
        dc -= 1;
    else
        dc += 1;
    return a + dr * 5 + dc;
}

bool PlayingAgainst(const State& state, God god) {
    for (size_t p = 0; p < state.gods.size(); p++)
        if (p != state.player && state.gods[p] == god) return true;
    return false;
}

bool AreNear(Coord a, Coord b) {
    int dc = int(a % 5) - int(b % 5);
    int dr = int(a / 5) - int(b / 5);
    return a != b && -1 <= dc && dc <= 1 && -1 <= dr && dr <= 1;
}

bool OnOppositeEdges(Coord a, Coord b) {
    int dc = int(a % 5) - int(b % 5);
    int dr = int(a / 5) - int(b / 5);
    return dc == -4 || dc == 4 || dr == -4 || dr == 4;
}

bool Enemy(const State& state, Cell d) {
    int pd = Player(d);
    return pd != -1 && pd != state.player;
}

// atomic actions:
// move worker (and optionally swap (apollo) or push (minotaur) or remove (bia))
// force worker (nemesis / charon / scylla / proteus)
// add worker (setup or later)
// remove worker (hydra, theseus)
// add dome
// add block (can build over opponent if medussa)
// remove block
// add token (also used for rotating wind token)

enum class ActionType : char {
    MoveWorker,
    ForceWorker,
    AddWorker, /*RemoveWorker,*/
    BuildDome,
    BuildBlock,
    RemoveBlock,
    AddToken
};
struct Action {
    ActionType type;
    char extra;  // worker/token when placing a new one, player that wins or looses
    Coord src, dest;
};

string Print(Action action) {
    // TODO
    return "";
}

optional<State> MoveWorker(const State& state, Coord src, Coord dest) {
    if (state.setup || src < 0 || src >= 25 || dest < 0 || dest >= 25 || !AreNear(src, dest)) return nullopt;

    Cell c = state[src];
    Cell d = state[dest];
    State s = state;
    God god = state.gods[state.player];

    if (Player(c) != state.player || Player(d) == state.player || Dome(d)) return nullopt;

    // Prometheus
    if (god != God::Prometheus && state.builds > 0) return nullopt;
    if (god == God::Prometheus && state.builds > 1) return nullopt;
    if (god == God::Prometheus && d.level > c.level && state.builds > 0) return nullopt;

    // Artemis
    if (god == God::Artemis) {
        if (state.moves > 1) return nullopt;
        if (state.moves > 0 && (state.lastMove != src || state.lastMoveSrc == dest)) return nullopt;
    }

    // Pegasus
    if (d.level - c.level > 1 && god != God::Pegasus) return nullopt;

    // Hades
    if (d.level < c.level && PlayingAgainst(state, God::Hades)) return nullopt;

    // Hermes
    if (god == God::Hermes) {
        if (d.level != c.level && !state.hermes_slow_move && state.moves > 0) return nullopt;
        if (state.hermes_slow_move) return nullopt;
        if (d.level != c.level) s.hermes_slow_move = true;
    }

    // Minotaur
    if (god == God::Minotaur) {
        int push = CellInDirection(src, dest);
        if (!Empty(d) && push != -1 && Empty(state[push])) swap(s[dest].figure, s[push].figure);
    }

    // Apollo
    if (!Empty(state[dest]) && god != God::Apollo) return nullopt;

    // Athena
    if (god == God::Athena && d.level > c.level) s.athenaMovedUp = true;
    if (d.level > c.level && state.athenaMovedUp && PlayingAgainst(state, God::Athena)) return nullopt;

    // Hera prevents victory when moving into perimeter space
    if (!PlayingAgainst(state, God::Hera) || !OnPerimeter(dest)) {
        if (d.level == 3 && c.level == 2) s.victory = true;

        // Pan
        if (god == God::Pan && d.level - c.level <= -2) s.victory = true;

        // Eros
        if (god == God::Eros && d.level == 1)
            for (Coord ie : CellsAround(dest))
                if (state[ie].level == 1 && Player(state[ie]) == state.player) s.victory = true;
    }

    // Bia
    if (god == God::Bia) {
        int target = CellInDirection(src, dest);
        if (target != -1) {
            Cell e = state[target];
            if (Enemy(state, e)) {
                s[target].figure = ' ';
                // TODO
                // if (!AnyWorker(s))
                //	s.victory = true;
            }
        }
    }

    s.lastMoveSrc = src;
    s.lastMove = dest;
    s.moves += 1;
    swap(s[src].figure, s[dest].figure);
    return s;
}

optional<State> ForceWorker(const State& state, Coord ia, Coord ib) {
    if (state.setup || ia < 0 || ia >= 25 || ib < 0 || ib >= 25) return nullopt;

    Cell a = state[ia];
    Cell b = state[ib];
    if (Dome(a) || Dome(b)) return nullopt;

    return nullopt;
    /*	// verify that ability allows forcing: Scylla, Proteus, Nemesis
            God god = state.gods[state.player];
            // swap workers on src and dest
            State s = state;
            swap(s[action.src].figure, s[action.dest].figure);
            return s;*/
}

optional<State> AddWorker(const State& state, Coord dest, char figure) {
    if (!state.setup || dest < 0 || dest >= 25 || !Empty(state[dest])) return nullopt;

    int male = 0, female = 0;
    Coord m = -1, f = -1;
    for (Coord i = 0; i < 25; i++) {
        char c = state[i].figure;
        if (c == 'A' + state.player) {
            male += 1;
            m = i;
        }
        if (c == 'a' + state.player) {
            female += 1;
            f = i;
        }
    }

    if (IsMale(figure) && male > 0) return nullopt;
    if (IsFemale(figure) && female > 0) return nullopt;
    if (Player(figure) != state.player) return nullopt;

    God god = state.gods[state.player];

    if (god == God::Eros) {
        Coord other = (m != -1) ? m : f;
        if (!OnOppositeEdges(m, dest)) return nullopt;
    }

    State s = state;
    s[dest].figure = figure;
    return s;
}

bool HasFiveCompleteTowers(const State& s) {
    int count = 0;
    for (int row = 0; row < 5; row++)
        for (int col = 0; col < 5; col++)
            if (auto c = s.cell[row][col]; Dome(c) && c.level == 3)
                if (++count == 5) return true;
    return false;
}

bool IsNearLimus(const State& state, Coord ie) {
    if (!PlayingAgainst(state, God::Limus)) return false;
    for (int ia : CellsAround(ie)) {
        int pa = Player(state[ia]);
        if (pa != -1 && pa != state.player && state.gods[pa] == God::Limus) return true;
    }
    return false;
}

optional<State> BuildDome(const State& state, Coord src, Coord dest) {
    if (state.setup || src < 0 || src >= 25 || dest < 0 || dest >= 25 || !AreNear(src, dest)) return nullopt;

    God god = state.gods[state.player];
    Cell s = state[src];
    Cell d = state[dest];

    if (!Empty(d) || Player(s) != state.player) return nullopt;

    State m = state;
    m.builds += 1;
    m.lastBuild = dest;
    m[dest].figure = ')';

    switch (god) {
        case God::Hermes:
            if (d.level < 3 || state.builds > 0) return nullopt;
            break;

        case God::Chronus:
            if (d.level < 3 || state.builds > 0 || state.lastMove != src) return nullopt;
            m.victory = HasFiveCompleteTowers(m);
            break;

        case God::Selene:
            if (state.lastBuild != -1) return nullopt;
            if (IsMale(s) && state.lastMove != src) return nullopt;
            break;

        case God::Atlas:
            if (state.builds > 0 || state.lastMove != src) return nullopt;
            break;

        case God::Demeter:
            if (d.level < 3 || state.builds >= 2 || state.lastMove != src) return nullopt;
            break;

        case God::Morpheus:
            if (d.level < 3 || state.morpheus_blocks == 0 || state.lastMove != src) return nullopt;
            m.morpheus_blocks -= 1;
            break;

        default:
            if (d.level < 3 || state.builds > 0 || state.lastMove != src) return nullopt;
    }

    if (d.level < 3 && IsNearLimus(state, dest)) return nullopt;
    if (PlayingAgainst(state, God::Chronus) && HasFiveCompleteTowers(m)) return nullopt;
    return m;
}

optional<State> BuildBlock(const State& state, Coord src, Coord dest) {
    if (state.setup || src < 0 || src >= 25 || dest < 0 || dest >= 25 || !AreNear(src, dest)) return nullopt;

    God god = state.gods[state.player];
    Cell s = state[src];
    Cell d = state[dest];

    if (!Empty(d) || Player(s) != state.player) return nullopt;

    if (IsNearLimus(state, dest)) return nullopt;

    State m = state;
    switch (god) {
        case God::Demeter:
            if (src != state.lastMove || state.builds > 1 || dest == state.lastBuild) return nullopt;
            break;

        case God::Hephaestus:
            if (src != state.lastMove || state.builds > 1 || (state.builds > 0 && dest != state.lastBuild))
                return nullopt;
            break;

        case God::Prometheus:
            if (state.moves == 0 && state.builds > 1) return nullopt;
            if (state.moves > 0 && src != state.lastMove && state.builds > 0) return nullopt;
            break;

        case God::Morpheus:
            if (src != state.lastMove || state.moves == 0 || state.morpheus_blocks == 0) return nullopt;
            m.morpheus_blocks -= 1;
            break;

        case God::Poseidon:
            break;

        case God::Hermes:
            if (state.builds > 0) return nullopt;
            break;

        default:
            if (src != state.lastMove || state.builds > 0) return nullopt;
    }

    m.builds += 1;
    m.lastBuild = dest;
    m[dest].level += 1;
    return m;
}

optional<State> RemoveBlock(const State& state, Coord dest) {
    if (dest < 0 || dest >= 25) return nullopt;

    God god = state.gods[state.player];
    if (god != God::Ares || state.builds == 0) return nullopt;

    Cell d = state[dest];
    if (!Empty(d) || d.level == 0) return nullopt;

    if (state.ares_removed_blocks > 0) return nullopt;

    // TODO check that there is unmoved worker of Ares near dest
    // TODO check that Ares can only remove one block per turn
    State m = state;
    m[dest].figure = ' ';
    m[dest].level -= 1;
    m.ares_removed_blocks += 1;
    return m;
}

optional<State> AddToken(const State& state, Coord dest, char token) { return nullopt; }

// returns nullopt if not allowed
optional<State> ExecuteAction(const State& state, Action action) {
    switch (action.type) {
        case ActionType::MoveWorker:
            return MoveWorker(state, action.src, action.dest);

        case ActionType::ForceWorker:
            return ForceWorker(state, action.src, action.dest);

        case ActionType::AddWorker:
            return AddWorker(state, action.dest, action.extra);

        case ActionType::BuildDome:
            return BuildDome(state, action.src, action.dest);

        case ActionType::BuildBlock:
            return BuildBlock(state, action.src, action.dest);

        case ActionType::RemoveBlock:
            return RemoveBlock(state, action.dest);

        case ActionType::AddToken:
            return AddToken(state, action.dest, action.extra);
    }
    return nullopt;
}

using Actions = vector<Action>;

#define ADD_STATE(s)         \
    if (s.victory) {         \
        states.clear();      \
        states.push_back(s); \
        return;              \
    } else                   \
        states.push_back(s);

void GenerateHermesFastMoves(const State& state, vector<State>& states) {
    // Find all horizontal possible moves for all builders at the same time
    static_vector<Coord, 25> builders;
    for (int row = 0; row < 5; row++)
        for (int col = 0; col < 5; col++)
            if (Player(state.cell[row][col]) == state.player) builders.push_back(row * 5 + col);

    if (builders.size() == 0 || builders.size() > 2) THROW(runtime_error);
    int sa = builders[0];
    int sb = (builders.size() == 2) ? builders[1] : 0;

    State base = state;
    char ca = ' ', cb = ' ';
    swap(base[sa].figure, ca);
    swap(base[sb].figure, cb);

    small_bfs<int> bfs(625);
    bfs.add(sa * 25 + sb, sa * 25 + sb);
    for (int e : bfs) {
        int ia = e / 25;
        int ib = e % 25;

        State s = base;
        s[ia].figure = ca;
        s[ib].figure = cb;
        s.lastMove = ia;
        ADD_STATE(s);

        // TODO move ia and add to bfs
        for (int ic : CellsAround(ia))
            if (ic != ib) {
            }

        if (builders.size() == 2) {
            // TODO move ib and add to bfs
        }
    }
}

bool InOrder(int a, int b, int c) { return a <= b && b <= c; }

constexpr int MAX_BUILDERS = 4;

bool AnyEnemy(const State& s) {
    for (int row = 0; row < 5; row++)
        for (int col = 0; col < 5; col++) {
            int p = Player(s.cell[row][col]);
            if (p != -1 && p != s.player) return true;
        }
    return false;
}

bool Wins(vector<State>& states) { return states.size() == 1 && states[0].victory; }

Coord Opposite(Coord a) { return (4 - a / 5) * 5 + (4 - a % 5); }

bool LessCells(const State& a, const State& b) {
    for (int i = 0; i < 25; i++)
        if (a[i] != b[i]) return a[i] < b[i];
    return false;
}

bool EqualCells(const State& a, const State& b) { return a.cell == b.cell; }

void Deduplicate(vector<State>& turns) { remove_dups(turns, LessCells, EqualCells); }

// TODO This is complicated. Could just use output vector of states as bfs queue with additional state hash set
void GenerateTritonMove(const State& state, int minJump, int maxJump, vector<State>& states) {
    for (int ic : PlayerCells(state)) {
        // TODO Triton can move back to initial cell if it jumps to perimeter first (it can also win this way)
        small_bfs<int> bfs(25);
        Cell c = state[ic];
        for (int id : CellsAround(ic)) {
            Cell d = state[id];

            // Important to check bfs.visited last as some level can only be reachable from certain direction
            if (!Empty(d) || !InOrder(minJump, d.level - c.level, maxJump) || bfs.visited[id]) continue;

            bfs.add(id, id);
            State s = state;
            // TODO check winning condition (and Hera)
            swap(s[ic].figure, s[id].figure);
            ADD_STATE(s);

            for (int ie : bfs)
                if (Cell e = state[ie]; OnPerimeter(ie))
                    for (int ib : CellsAround(ie))
                        if (Cell b = state[ib]; Empty(b) && InOrder(minJump, b.level - e.level, maxJump)) {
                            bfs.add(ib, ib);
                            s = state;
                            // TODO check winning condition (and Hera)
                            swap(s[ic].figure, s[ib].figure);
                            ADD_STATE(s);
                        }
        }
    }
}

// TODO This is complicated, and wrong! Artemis might take two moves to move up one cell from current cell.
void GenerateArtemisMoves(const State& state, int minJump, int maxJump, vector<State>& states) {
    // perform one or two moves with a single builder, but not back to initial cell
    for (auto ic : PlayerCells(state)) {
        Cell c = state[ic];

        for (auto id : CellsAround(ic)) {
            Cell d = state[id];
            if (!Empty(d) || !InOrder(minJump, d.level - c.level, maxJump)) continue;

            State s = MoveBuilder(state, ic, id);
            ADD_STATE(s);

            // second move (in the same direction)
            auto ie = CellInDirection(ic, id);
            if (ie == -1) continue;
            Cell e = state[ie];
            if (!Empty(e) || !InOrder(minJump, e.level - d.level, maxJump)) continue;

            State u = MoveBuilder(s, id, ie);
            ADD_STATE(u);
        }

        for (auto id : CellsAroundL2(ic)) {
            Cell d = state[id];
            if (!Empty(d)) continue;

            int dr = id / 5 - ic / 5;
            int dc = id % 5 - ic % 5;

            Coord ia = LMiddle1(ic, dr, dc);
            Cell a = state[ia];
            if (Empty(a) && InOrder(minJump, a.level - c.level, maxJump) &&
                InOrder(minJump, d.level - a.level, maxJump)) {
                State s = MoveBuilder(state, ic, ia);
                s = MoveBuilder(s, ia, id);
                ADD_STATE(s);
                continue;
            }

            Coord ib = LMiddle2(ic, dr, dc);
            Cell b = state[ib];
            if (Empty(b) && InOrder(minJump, b.level - c.level, maxJump) &&
                InOrder(minJump, d.level - b.level, maxJump)) {
                State s = MoveBuilder(state, ic, ib);
                s = MoveBuilder(s, ib, id);
                ADD_STATE(s);
            }
        }
    }
}

void GenerateOneMove(const State& state, vector<State>& states) {
    God god = state.gods[state.player];

    int minJump = -3;
    if (god != God::Hades && contains(state.gods, God::Hades)) minJump = 0;

    int maxJump = 1;
    if (god == God::Pegasus) maxJump = 3;
    if (god == God::Prometheus && state.lastBuild != -1) maxJump = 0;
    if (state.athenaMovedUp && god != God::Athena) maxJump = 0;

    if (god == God::Hermes) GenerateHermesFastMoves(state, states);

    if (god == God::Triton) {
        GenerateTritonMove(state, minJump, maxJump, states);
        return;
    }

    if (god == God::Artemis) {
        GenerateArtemisMoves(state, minJump, maxJump, states);
        return;
    }

    for (int ic : PlayerCells(state))
        for (int id : CellsAround(ic)) {
            Cell c = state[ic];
            Cell d = state[id];

            if (Player(d) == state.player || Dome(d) || !InOrder(minJump, d.level - c.level, maxJump)) continue;

            // horizontal case is covered separately
            if (god == God::Hermes && d.level == c.level) continue;

            int ie = -1;
            if (god == God::Minotaur) {
                ie = CellInDirection(ic, id);
                if (!Empty(d) && (ie == -1 || !Empty(state[ie]))) continue;
            } else if (god != God::Apollo && !Empty(d))
                continue;

            State s = state;
            s.lastMove = id;
            if (god == God::Athena && d.level > c.level) s.athenaMovedUp = true;
            if (god == God::Pan && d.level - c.level <= -2) s.victory = true;
            if (d.level == 3 && c.level == 2) {
                s.victory = true;
                // Hera prevents others from winning by climbing on a perimeter tower
                if (god != God::Hera && contains(state.gods, God::Hera) && !OnPerimeter(id)) s.victory = false;
            }
            if (god == God::Eros && d.level == 1)
                for (Coord ie : CellsAround(id))
                    if (state[ie].level == 1 && Player(state[ie]) == state.player) {
                        s.victory = true;
                        // Hera prevents Eros's winning condition by moving to perimeter
                        if (contains(state.gods, God::Hera) && !OnPerimeter(id)) s.victory = false;
                    }

            if (god == God::Minotaur) swap(s[ie].figure, s[id].figure);
            swap(s[ic].figure, s[id].figure);
            ADD_STATE(s);
        }
}

// build in one cell (hephaestus can build twice in one cell)
void GenerateOneBuild(int builder, const State& state, vector<State>& states, bool seleneDome = false) {
    God god = state.gods[state.player];
    bool containsChronus = contains(state.gods, God::Chronus);
    bool containsLimus = contains(state.gods, God::Limus);

    for (int ie : CellsAround(builder, god == God::Zeus)) {
        Cell e = state[ie];
        if (!Empty(e) || (state.lastBuild == ie && god == God::Demeter)) continue;

        if (e.level == 3) {
            if (seleneDome) continue;
            // build complete tower
            State s = state;
            s[ie].figure = ')';
            s.lastBuild = ie;

            if (containsChronus && HasFiveCompleteTowers(s)) {
                if (god == God::Chronus)
                    s.victory = true;
                else
                    continue;
            }

            ADD_STATE(s);
            continue;
        }

        if (containsLimus && god != God::Limus && IsNearLimus(state, ie)) continue;

        if (god == God::Atlas || seleneDome) {
            // build early dome
            State s = state;
            s[ie].figure = ')';
            s.lastBuild = ie;
            ADD_STATE(s);
        }

        if (!seleneDome) {
            // build floor
            State s = state;
            s[ie].level += 1;
            s.lastBuild = ie;
            ADD_STATE(s);

            if (e.level < 2 && god == God::Hephaestus) {
                s[ie].level += 1;
                ADD_STATE(s);
            }
        }
    }
}

void GenerateTurns(State state, vector<State>& turns) {
    God god = state.gods[state.player];
    if (state.setup) {
        // TODO if state.player == 0 then use a smaller set of initial positions (ie. symmetry)
        if (god == God::Eros) {
            for (Coord ia = 0; ia < 25; ia++)
                if (OnPerimeter(ia) && Empty(state[ia])) {
                    Coord ib = Opposite(ia);
                    if (ia < ib && Empty(state[ib])) {
                        State s = state;
                        s[ia].figure = 'A' + state.player;
                        s[ib].figure = 'a' + state.player;
                        s.player += 1;
                        if (s.player == s.gods.size()) {
                            s.player = 0;
                            s.setup = false;
                        }
                        turns.push_back(s);
                    }
                }
            return;
        }

        for (Coord ia = 0; ia < 25; ia++)
            if (Empty(state[ia]))
                for (Coord ib = ia + 1; ib < 25; ib++)
                    if (Empty(state[ib])) {
                        State s = state;
                        s[ia].figure = 'A' + state.player;
                        s[ib].figure = 'a' + state.player;
                        s.player += 1;
                        if (s.player == s.gods.size()) {
                            s.player = 0;
                            s.setup = false;
                        }
                        turns.push_back(s);

                        if (god == God::Selene) {
                            swap(s[ia].figure, s[ib].figure);
                            turns.push_back(s);
                        }
                    }
        return;
    }

    if (god == God::Athena) state.athenaMovedUp = false;
    state.lastMove = -1;
    state.lastBuild = -1;
    turns.push_back(state);
    thread_local vector<State> temp;

    // perform optional build if prometheus
    if (god == God::Prometheus)
        for (auto ia : PlayerCells(state)) GenerateOneBuild(ia, state, turns);

    // perform one move
    temp.clear();
    swap(temp, turns);
    for (const State& m : temp) GenerateOneMove(m, turns);
    if ((turns.empty() && god != God::Hermes) || Wins(turns)) return;

    // perform one build
    temp.clear();
    swap(temp, turns);
    for (const State& m : temp)
        if (god == God::Hermes)
            for (Coord builder : PlayerCells(state)) GenerateOneBuild(builder, m, turns);
        else
            GenerateOneBuild(m.lastMove, m, turns);

    if (god == God::Selene)
        // Female builder can build a dome at any level even if it didn't move (instead of normal build)
        for (const State& m : temp)
            for (Coord builder : PlayerCells(state))
                if (IsFemale(m[builder])) GenerateOneBuild(builder, m, turns, /*seleneDome=*/true);

    if (turns.empty() || Wins(turns)) return;

    // perform second optional build if demeter
    if (god == God::Demeter) {
        temp.clear();
        swap(temp, turns);
        for (const State& m : temp) {
            turns.push_back(m);  // second build is optional
            GenerateOneBuild(m.lastMove, m, turns, m.lastBuild);
        }
        Deduplicate(turns);
    }

    // perform three additional builds with unmoved builder on the ground level if poseidon
    if (god == God::Poseidon) {
        for (int i = 0; i < 3; i++) {
            temp.clear();
            swap(temp, turns);
            for (const State& m : temp) {
                turns.push_back(m);  // building is optional
                for (auto ia : PlayerCells(m))
                    if (ia != m.lastMove && m[ia].level == 0) GenerateOneBuild(ia, m, turns);
            }
            Deduplicate(turns);
        }
    }

    if (god == God::Medusa)
        for (State& m : turns)
            for (auto ia : PlayerCells(m))
                for (auto ib : CellsAround(ia)) {
                    Cell a = m[ia];
                    Cell& b = m[ib];
                    if (a.level > b.level && Enemy(m, b)) {
                        b.level += 1;
                        b.figure = ' ';
                        m.lastBuild = ib;
                        if (!AnyEnemy(m)) {
                            m.victory = true;
                            turns[0] = m;
                            turns.resize(1);
                            return;
                        }
                    }
                }

    // optionally destroy a non-dome unoccupied level near unmoved worker if Ares
    if (god == God::Ares) {
        temp.clear();
        swap(temp, turns);
        for (const State& m : temp) {
            turns.push_back(m);
            for (auto ia : PlayerCells(m))
                if (ia != m.lastMove)
                    for (auto ib : CellsAround(ia))
                        if (m[ib].level > 0 && Empty(m[ib])) {
                            State s = m;
                            s[ib].level -= 1;
                            turns.push_back(s);
                        }
        }
    }

    // TODO assert no duplicates

    // end turn
    int next = NextPlayer(state);
    for (State& s : turns) s.player = next;
}

optional<State> EndTurn(const State& state) {
    God god = state.gods[state.player];
    State m = state;

    if (state.setup) {
        // TODO verify both workers have been placed
    }

    // check if player completed its turn
    if (state.moves == 0 && god != God::Hermes) return nullopt;
    if (state.builds == 0 && god != God::Morpheus) return nullopt;

    // auto-stone enemies around Medussa on lower level
    if (god == God::Medusa)
        for (auto ia : PlayerCells(m))
            for (auto ib : CellsAround(ia)) {
                Cell a = m[ia];
                Cell& b = m[ib];
                if (a.level > b.level && Enemy(m, b)) {
                    b.level += 1;
                    b.figure = ' ';
                    if (!AnyEnemy(m)) {
                        m.victory = true;
                        break;
                    }
                }
            }

    // switch state for next player
    m.player = NextPlayer(state);
    if (m.setup && m.player == 0) m.setup = false;
    m.lastMoveSrc = -1;
    m.lastMove = -1;
    m.lastBuild = -1;
    m.builds = 0;
    m.moves = 0;
    m.hermes_slow_move = false;
    m.ares_removed_blocks = 0;
    god = m.gods[m.player];

    // if next player is Morpheus give him a block
    if (god == God::Morpheus) m.morpheus_blocks += 1;

    if (god == God::Athena) m.athenaMovedUp = false;
    return m;
}

void GenerateTurns(const State& state, vector<pair<State, Actions>>& states) {
    God god = state.gods[state.player];

    // Move or build abilities
    if (god == God::Artemis) {
        // TODO add all possible moves from state to states
        // TODO for each state in states, generate all possible second moves (ignoring dups)
        // TODO end turn on every state
        return;
    }

    if (god == God::Hermes) {
        // TODO add state to states
        // TODO for each state in states (BFS), generate all possible moves (ignoring dups)

        // TODO for each s in states generate all possible builds (not-optional)
        // TODO end turn on every state
        return;
    }

    if (god == God::Prometheus) {
        // TODO add all possible moves and builds from state to states
        // TODO move after build (not optional)
        // TODO build after move (not optional)
        // TODO end turn on every state
        return;
    }

    if (god == God::Morpheus) {
        // TODO add all possible moves from state to states
        // TODO for each state in states (BFS), generate all possible builds (ignoring dups)
        // TODO end turn on every state
        return;
    }

    if (god == God::Selene) {
        // TODO end turn on every state
        return;
    }

    if (god == God::Triton) {
        // TODO add all possible moves from state to states
        // TODO for each state in states (BFS), if on perimeter, generate all possible moves (ignoring dups)

        // TODO for each s in states generate all possible builds
        // TODO end turn on every state
        return;
    }

    if (god == God::Zeus) {
        // TODO end turn on every state
        return;
    }

    // everyone else
    // TODO generate all possible moves from state
    // TODO for each s in states: generate all possible builds

    // End of turn abilities
    if (god == God::Demeter) {
        // TODO for each s in states: generate optional second build on a different space
    }

    if (god == God::Ares) {
        // TODO for each s in states: optionally remove one unoccupied block
    }

    if (god == God::Hephaestus) {
        // TODO for each s in states: generate optional second build on the same space
    }

    if (god == God::Poseidon) {
        // TODO 3 optional builds with unmoved worker on ground level (ignoring dups)
    }

    // TODO end turn on every state
}

// TODO print possible moves side by side to save space
// TODO print last build and last move in different color
void Print(const State& state) {
    print("player %s", int(state.player));
    if (state.golden_fleece != God::None) {
        print(" golden fleece %s", magic_enum::enum_name(state.golden_fleece));
    } else {
        for (God a : state.gods) print(" %s", magic_enum::enum_name(a));
    }
    if (state.wind != -1) print(", wind %d", state.wind);
    if (state.arrow != -1) print(", arrow %d", state.arrow);
    if (state.athenaMovedUp) print(", athena moved up");
    if (state.victory) print(", victory!");
    print("\n");

    for (auto& row : state.cell) {
        for (auto& c : row) {
            char code[] = {' ', ' ', ' ', ' ', ' ', '\0'};
            char* p = code;
            for (int i = 0; i < c.level; i++) *p++ = ']';

            int tokens = 0;
            if (c.coin) {
                *p = '$';
                tokens += 1;
            }
            if (c.whirlpool) {
                *p = '*';
                tokens += 1;
            }
            if (c.abyss) {
                *p = 'O';
                tokens += 1;
            }
            if (c.fate) {
                *p = '#';
                tokens += 1;
            }
            if (tokens == 1) p++;
            if (tokens > 1) *p++ = '?';

            *p++ = c.figure;
            if (code[0] == ' ') code[0] = '.';
            print("%s", code);
        }
        print("\n");
    }
    print("\n");
}

// TODO add option to reset to start of turn
optional<State> Human(const State& state) {
    vector<State> moves;
    GenerateTurns(state, moves);
    print("total posible moves %s\n", moves.size());
    if (moves.size() == 0) return nullopt;
    for (int i = 0; i < moves.size(); i++) {
        print("possible move %s\n", i);
        Print(moves[i]);
    }
    while (true) {
        print("choose your move > ");
        int d = -1;
        if (scanf("%d", &d) == 1 && 0 <= d && d < moves.size()) return moves[d];
    }
}

optional<State> RandBot(const State& state) {
    vector<State> moves;
    GenerateTurns(state, moves);
    if (moves.size() == 0) return nullopt;
    return moves[rand() % moves.size()];
}

template <typename T, typename E>
T bestElement(cspan<T> data, const E& scoreFunc) {
    vector<size_t> best = {0};
    auto bestScore = scoreFunc(data[0]);
    for (size_t i = 1; i < data.size(); i++) {
        auto score = scoreFunc(data[i]);
        if (score > bestScore) {
            bestScore = score;
            best = {i};
        } else if (score == bestScore)
            best.push_back(i);
    }
    return data[best[rand() % best.size()]];
}

optional<State> GreedyBot(const State& state) {
    vector<State> moves;
    GenerateTurns(state, moves);
    if (moves.size() == 0) return nullopt;
    if (moves.size() == 1 && moves[0].victory) return moves[0];

    return bestElement(cspan<State>(moves), [&state](const State& move) {
        double score = 0;
        for (int ic : PrevPlayerCells(move)) {
            Cell c = move[ic];
            score += 100 * c.level * c.level;
            for (int id : CellsAround(ic)) {
                Cell d = move[id];
                if (Empty(d)) score += d.level * d.level;
                // blocking with domes
                Cell d2 = state[id];
                if (Dome(d) && !Dome(d2) && d.level == 3) {
                    bool blocking = false;
                    for (int ie : CellsAround(id)) {
                        Cell e = move[ie];
                        if (!Empty(e) && e.level == 2 && Player(e) != state.player) blocking = true;
                    }
                    if (blocking) score += 100000;
                }
            }
        }
        for (int ic : PlayerCells(move)) {
            Cell c = move[ic];
            score -= 100 * c.level * c.level;
            for (int id : CellsAround(ic)) {
                Cell d = move[id];
                if (Empty(d)) {
                    if (d.level == 3) score = -INF;
                    score -= d.level * d.level;
                }
            }
        }
        return score;
    });
}

double Heuristic(const State& state, const State& move) {
    double score = 0;
    for (int ic : PrevPlayerCells(move)) {
        Cell c = move[ic];
        score += 100 * c.level * c.level;
        for (int id : CellsAround(ic)) {
            Cell d = move[id];
            if (Empty(d)) score += d.level * d.level;
            // blocking with domes
            Cell d2 = state[id];
            if (Dome(d) && !Dome(d2) && d.level == 3) {
                bool blocking = false;
                for (int ie : CellsAround(id)) {
                    Cell e = move[ie];
                    if (!Empty(e) && e.level == 2 && Player(e) != state.player) blocking = true;
                }
                if (blocking) score += 100000;
            }
        }
    }
    for (int ic : PlayerCells(move)) {
        Cell c = move[ic];
        score -= 100 * c.level * c.level;
        for (int id : CellsAround(ic)) {
            Cell d = move[id];
            if (Empty(d)) {
                if (d.level == 3) score = -INF;
                score -= d.level * d.level;
            }
        }
    }
    return score;
}

// TODO generalize for multiple state updates from a single player
double SubTreeScore(const State& prev, const State& state, int depth, bool maxi) {
    if (depth == 0) return Heuristic(prev, state) * (maxi ? 1 : -1);

    vector<State> moves;
    GenerateTurns(state, moves);
    if (moves.size() == 0) return maxi ? -INF : INF;
    if (moves.size() == 1 && moves[0].victory) return maxi ? INF : -INF;

    double best = (depth % 2) ? INF : -INF;
    for (const State& m : moves) {
        double s = SubTreeScore(state, m, depth - 1, !maxi);
        best = maxi ? max(best, s) : min(best, s);
    }
    return best;
}

// TODO generalize for multiple state updates from a single player
optional<State> BruteBot(const State& state) {
    vector<State> moves;
    GenerateTurns(state, moves);
    if (moves.size() == 0) return nullopt;
    if (moves.size() == 1 && moves[0].victory) return moves[0];

    return bestElement(cspan<State>(moves), [&state](const State& m) { return SubTreeScore(state, m, 2, false); });
}

using Strategy = std::function<optional<State>(const State&)>;

double RelativeSkill2(int games, vector<God> gods, Strategy a, Strategy b) {
    vector<Strategy> strategy = {a, b};
    vector<int> wins;
    wins.resize(gods.size(), 0);
    for (int i = 0; i < games; i++) {
        State state;
        for (int j = 0; j < state.gods.size(); j++) state.gods[j] = (j < gods.size()) ? gods[j] : God::Dead;

        while (true) {
            auto s = strategy[state.player](state);
            if (!s.has_value()) {
                // TODO 2 player assumption
                wins[1 - state.player] += 1;
                break;
            }
            state = *s;
            if (state.victory) {
                wins[state.player] += 1;
                break;
            }
        }
    }
    return double(wins[0]) / games;
}

// TODO make it multi-threaded!
double RelativeSkill(int games, vector<God> gods, Strategy a, Strategy b) {
    return (RelativeSkill2(games / 2, gods, a, b) + 1 - RelativeSkill2(games / 2, {gods[1], gods[0]}, b, a)) / 2;
}

void SingleMatch(cspan<God> gods, Strategy a, Strategy b, bool verbose = true) {
    State state;
    for (int j = 0; j < state.gods.size(); j++) state.gods[j] = (j < gods.size()) ? gods[j] : God::Dead;
    if (verbose) Print(state);

    vector<Strategy> strategy = {a, b};
    while (!state.victory) {
        auto s = strategy[state.player](state);
        if (!s.has_value()) break;
        state = *s;
        if (verbose) Print(state);
    }
}

TEST_CASE("Santorini Apollo", "[santorini]") {
    print("running\n");
    State state =
        ".b. . . . "
        ". 2 2 2 . "
        ". 2 .a2B. "
        ". 2 2A1). "
        ". . 2)2 . "sv;
    state.gods = {God::Apollo, God::None};
    state.setup = false;
    Print(state);

    vector<State> moves;
    GenerateTurns(state, moves);
    // for (const State& m : moves)
    //	Print(m);
}

inline God operator++(God& x) { return x = God(int(x) + 1); }

TEST_CASE("All pairs matchups", "[santorini]") {
    for (God a : magic_enum::enum_values<God>())
        for (God b : magic_enum::enum_values<God>())
            if (a != God::Dead && b != God::Dead) {
                print("%s vs %s\n", magic_enum::enum_name(a), magic_enum::enum_name(b));
                SingleMatch({a, b}, RandBot, RandBot, false);
                SingleMatch({b, a}, RandBot, RandBot, false);
                SingleMatch({a, b}, GreedyBot, GreedyBot, false);
                SingleMatch({b, a}, GreedyBot, GreedyBot, false);
            }
}

int main(int argc, char* argv[]) {
    InitSegvHandler();
    srand(time(0));

    if (argc == 2 && strcmp("-t", argv[1]) == 0) {
        char* test_argv[1] = {argv[0]};
        return Catch::Session().run(1, test_argv);
    }

    SingleMatch({God::Ares, God::None}, Human, RandBot);

    vector<God> gods = {God::None, God::None};
    print("Greedy - Rand %s\n", RelativeSkill(10000, gods, GreedyBot, RandBot));
    print("Brute - Rand %s\n", RelativeSkill(10, gods, BruteBot, RandBot));
    print("Brute - Greedy %s\n", RelativeSkill(10, gods, BruteBot, GreedyBot));
    return 0;
}
