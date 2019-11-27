#include <core/format.h>
#include <core/each.h>
#include <core/dynamic_array.h>
#include <core/exception.h>
#include <core/util.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum class God : char { None, Dead,
	Apollo, Artemis, Athena, Atlas, Demeter, Hephaestus, /*partial*/Hermes, Minotaur, Pan, Prometheus,
	Chronus, Hera, Limus, Medusa, Poseidon, Zeus };

struct Cell {
	bool dome = false;
	char tower = 0; // 0-3 height of tower
	char builder = ' '; // Aa - first player, Bb - second player, ... (uppercase male, lowercase female)
};

bool operator<(Cell a, Cell b) {
	if (a.dome != b.dome)
		return !a.dome && b.dome;
	if (a.tower != b.tower)
		return a.tower < b.tower;
	return a.builder < b.builder;
}

bool operator==(Cell a, Cell b) {
	return a.dome == b.dome && a.tower == b.tower && a.builder == b.builder;
}

bool operator!=(Cell a, Cell b) { return !(a == b); }

struct State {
	array<array<Cell, 5>, 5> cell;
	array<God, 2> gods;
	bool athenaMovedUp = false;
	char player = 0; // player about to play
	char lastMove = -1; // coordinate of builder who moved last
	char lastBuild = -1; // coordinate of last build
	bool victory = false;

	Cell operator[](int a) const { return cell[a / 5][a % 5]; }
	Cell& operator[](int a) { return cell[a / 5][a % 5]; }
};

auto Cells() {
	array<int, 25> data;
	for (int i = 0; i < data.size(); i++)
		data[i] = i;
	return data;
}

int Player(Cell cell) {
	auto builder = cell.builder;
	if ('a' <= builder && builder <= 'd')
		return builder - 'a';
	if ('A' <= builder && builder <= 'D')
		return builder - 'A';
	return -1;
}

int PrevPlayer(const State& state) {
	for (int i = 1; i < state.gods.size(); i++) {
		int p = (state.player + state.gods.size() - i) % state.gods.size();
		if (state.gods[p] != God::Dead)
			return p;
	}
	THROW(runtime_error);
}

int NextPlayer(const State& state) {
	for (int i = 1; i < state.gods.size(); i++) {
		int p = (state.player + i) % state.gods.size();
		if (state.gods[p] != God::Dead)
			return p;
	}
	THROW(runtime_error);
}

auto PrevPlayerCells(const State& state) {
	int p = PrevPlayer(state);
	static_vector<int, 4> data;
	for (int i = 0; i < 25; i++)
		if (Player(state[i]) == p)
			data.push_back(i);
	return data;
}

auto PlayerCells(const State& state) {
	static_vector<int, 4> data;
	for (int i = 0; i < 25; i++)
		if (Player(state[i]) == state.player)
			data.push_back(i);
	return data;
}

// TODO precompute
auto CellsAround(int a, bool includeCenter = false) {
	static_vector<int, 8> data;
	int row = a / 5;
	int col = a % 5;
	for (int dr = -1; dr <= 1; dr++)
		for (int dc = -1; dc <= 1; dc++) {
			if (dr == 0 && dc == 0 && !includeCenter)
				continue;
			int r = row + dr;
			int c = col + dc;
			if (0 <= r && r < 5 && 0 <= c && c < 5)
				data.push_back(r * 5 + c);
		}
	return data;
}

int CellInDirection(int a, int b) {
	int row = b / 5 + b / 5 - a / 5;
	int col = b % 5 + b % 5 - a % 5;
	return (0 <= row && row < 5 && 0 <= col && col < 5) ? row * 5 + col : -1;
}

// TODO precompute
bool OnPerimeter(int a) {
	int row = a / 5;
	int col = a % 5;
	return row == 0 || row == 4 || col == 0 || col == 4;
}

void GenerateOneMove(const State& state, vector<State>& moves, bool allowMoveUp, const State* forbidden) {
	God god = state.gods[state.player];
	if (god == God::Hermes) {
		// TODO find all possible moves for all builders
	}

	for (int ic : PlayerCells(state))
		for (int id : CellsAround(ic)) {
			Cell c = state[ic];
			Cell d = state[id];

			if (Player(d) == state.player)
				continue;
			if (d.dome)
				continue;

			// horizontal case is covered separately
			if (god == God::Hermes && d.tower == c.tower)
				continue;

			if (!allowMoveUp && d.tower > c.tower)
				continue;

			int maxJump = 1;
			// TODO move this into allowMoveUp
			if (state.athenaMovedUp && god != God::Athena)
				maxJump = 0;
			if (d.tower - c.tower > maxJump)
				continue;

			int ie = -1;
			if (god == God::Minotaur) {
				ie = CellInDirection(ic, id);
				if (d.builder != ' ' && (ie == -1 || state[ie].builder != ' '))
					continue;
			} else if (god != God::Apollo && d.builder != ' ')
				continue;

			State s = state;
			s.lastMove = id;
			if (d.tower - c.tower > 0 && god == God::Athena)
				s.athenaMovedUp = true;
			if (god == God::Pan && d.tower - c.tower <= -2)
				s.victory = true;
			if (c.tower < 3 && d.tower == 3)
				if (!contains(state.gods, God::Hera) || !OnPerimeter(id))
					s.victory = true;

			if (god == God::Minotaur)
				swap(s[ie].builder, s[id].builder);
			swap(s[ic].builder, s[id].builder);

			if (forbidden && s.cell == forbidden->cell)
				continue;
			moves.push_back(s);
		}
}

constexpr int MAX_BUILDERS = 4;

int CompleteTowerCount(const State& s) {
	int count = 0;
	for (int row = 0; row < 5; row++)
		for (int col = 0; col < 5; col++)
			if (s.cell[row][col].dome && s.cell[row][col].tower == 3)
				count += 1;
	return count;
}

int BuilderCount(const State& s, int player) {
	int count = 0;
	for (int row = 0; row < 5; row++)
		for (int col = 0; col < 5; col++)
			if (player == Player(s.cell[row][col]))
				count += 1;
	return count;
}

// build in one cell (hephaestus can build twice in one cell)
void GenerateOneBuild(int builder, const State& state, vector<State>& builds, int forbidden = -1) {
	God god = state.gods[state.player];
	for (int ie : CellsAround(builder, god == God::Zeus)) {
		Cell e = state[ie];
		if (e.dome || e.builder != ' ' || ie == forbidden)
			continue;

		bool nearLimus = false;
	    if (god != God::Limus && contains(state.gods, God::Limus))
			for (int ia : CellsAround(ie)) {
				int pa = Player(state[ia]);
				if (pa != -1 && state.gods[pa] == God::Limus) {
					nearLimus = true;
					break;
				}
			}

		if (e.tower == 3 || (god == God::Atlas && !nearLimus)) {
			State s = state;
			s[ie].dome = true;
			s.lastBuild = ie;

			if (e.tower == 3 && contains(state.gods, God::Chronus) && CompleteTowerCount(s) == 5) {
				for (size_t i = 0; i < state.gods.size(); i++)
					if (state.gods[i] == God::Chronus)
						s.player = i;
				s.victory = true;
			}
			builds.push_back(s);
		}

		if (e.tower < 3 && !nearLimus) {
			State s = state;
			s[ie].tower += 1;
			s.lastBuild = ie;
			builds.push_back(s);

			if (e.tower < 2 && god == God::Hephaestus) {
				s[ie].tower += 1;
				builds.push_back(s);
			}
		}
	}
}

bool FilterVictory(vector<State>& states) {
	for (State& s : states)
		if (s.victory) {
			states = { s };
			return true;
		}
	return false;
}

// TODO: initial figure placement should also be a turn!
void GenerateTurns(State state, vector<State>& turns) {
	God god = state.gods[state.player];
	if (god == God::Athena)
		state.athenaMovedUp = false;
	state.lastMove = -1;
	state.lastBuild = -1;
	turns.push_back(state);

	// perform optional build if prometheus
	if (god == God::Prometheus)
		for (int ia : PlayerCells(state))
			GenerateOneBuild(ia, state, turns);

	// perform one move
	vector<State> temp;
	swap(temp, turns);
	for (const State& m : temp)
		GenerateOneMove(m, turns, god != God::Prometheus || m.lastBuild == -1, nullptr);
	if (turns.size() == 0 && god != God::Hermes)
		return;
	if (FilterVictory(turns))
		return;

	// perform second optional move if artemis
	if (god == God::Artemis) {
		temp.clear();
	    swap(temp, turns);
		for (const State& m : temp)
			GenerateOneMove(m, turns, true, &state);
		if (FilterVictory(turns))
			return;
	}

	// perform one build
	temp.clear();
	swap(temp, turns);
	for (const State& m : temp)
		if (god == God::Hermes)
			for (int builder : PlayerCells(state))
				GenerateOneBuild(builder, m, turns);
		else
			GenerateOneBuild(m.lastMove, m, turns);
	if (turns.size() == 0)
		return;
	if (FilterVictory(turns))
		return;

	// perform second optional build if demeter
	if (god == God::Demeter) {
		temp.clear();
		swap(temp, turns);
		for (const State& m : temp) {
			turns.push_back(m); // second build is optional
			GenerateOneBuild(m.lastMove, m, turns, m.lastBuild);
		}
		if (FilterVictory(turns))
			return;
	}

	// perform three additional builds with unmoved builder on the ground level if poseidon
	if (god == God::Poseidon) {
		std::function<bool(const State&, const State&)> less = [](const State& a, const State& b) {
			for (int i = 0; i < 25; i++)
				if (a[i] != b[i])
					return a[i] < b[i];
			return false;
		};
		std::function<bool(const State&, const State&)> equal = [](const State& a, const State& b) {
			return a.cell == b.cell;
		};
		for (int i = 0; i < 3; i++) {
			temp.clear();
			swap(temp, turns);
			for (const State& m : temp) {
				turns.push_back(m); // building is optional
				for (int ia : PlayerCells(state))
					if (ia != m.lastMove && m[ia].tower == 0)
						GenerateOneBuild(ia, m, turns, -1);
				remove_dups(turns, less, equal);
			}
			if (FilterVictory(turns))
				return;
		}
	}

	// perform medusa's build
	if (god == God::Medusa) {
		for (State& m : turns) {
			for (int ia : PlayerCells(m))
				for (int ib : CellsAround(ia)) {
					Cell a = m[ia];
					Cell& b = m[ib];
					int pb = Player(b);
					if (a.tower > b.tower && pb != -1 && pb != state.player) {
						b.tower += 1;
						b.builder = ' ';
						m.lastBuild = ib;
						if (BuilderCount(m, pb) == 0) {
							// TODO 2 player assumption
							m.victory = true;
						}
					}
				}
		}
		if (FilterVictory(turns))
			return;
	}

	// end turn
	int next = NextPlayer(state);
	for (State& s : turns)
		s.player = next;
}

void Print(const State& state) {
	print("player %s", int(state.player));
	if (state.lastMove != -1)
		print(", move (%s %s)", state.lastMove / 5, state.lastMove % 5);
	if (state.lastBuild != -1)
		print(", build (%s %s)", state.lastBuild / 5, state.lastBuild % 5);
	if (state.athenaMovedUp)
		print(", athena moved up");
	if (state.victory)
		print(", victory!");
	print("\n");

	for (auto& row : state.cell) {
		for (auto& c : row) {
			char code[] = {' ', ' ', ' ', ' ', ' ', '\0'};
			char* p = code;
			for (int i = 0; i < c.tower; i++)
				*p++ = ']';
			if (c.dome)
				*p++ = ')';
			*p++ = c.builder;
			if (code[0] == ' ')
				code[0] = '.';
			print("%s", code);
		}
		print("\n");
	}
	print("\n");
}

void PlaceBuilder(State& state, char builder) {
	while (true) {
		Cell& c = state.cell[rand() % 5][rand() % 5];
		if (c.builder == ' ') {
			c.builder = builder;
			return;
		}
	}
}

// TODO print possible moves side by side to save space
optional<State> Human(const State& state) {
	vector<State> moves;
	GenerateTurns(state, moves);
	print("total posible moves %s\n", moves.size());
	if (moves.size() == 0)
		return nullopt;
	for (int i = 0; i < moves.size(); i++) {
		print("possible move %s\n", i);
		Print(moves[i]);
	}
	while (true) {
		print("choose your move > ");
		int d = -1;
		if (scanf("%d", &d) == 1 && 0 <= d && d < moves.size())
			return moves[d];
	}
}

optional<State> RandBot(const State& state) {
	vector<State> moves;
	GenerateTurns(state, moves);
	if (moves.size() == 0)
		return nullopt;
	return moves[rand() % moves.size()];
}

template<typename T, typename E>
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
	if (moves.size() == 0)
		return nullopt;
	if (moves.size() == 1 && moves[0].victory)
		return moves[0];

	return bestElement(cspan<State>(moves), [&state](const State& move){
		double score = 0;
		for (int ic : PrevPlayerCells(move)) {
			Cell c = move[ic];
			score += 100 * c.tower * c.tower;
			for (int id : CellsAround(ic)) {
				Cell d = move[id];
				if (d.builder == ' ' && !d.dome)
					score += d.tower * d.tower;
				// blocking with domes
				Cell d2 = state[id];
				if (d.dome && !d2.dome && d.tower == 3) {
					bool blocking = false;
					for (int ie : CellsAround(id)) {
						Cell e = move[ie];
						if (e.builder != ' ' && e.tower == 2 && Player(e) != state.player)
							blocking = true;
					}
					if (blocking)
						score += 100000;
				}
			}
		}
		for (int ic : PlayerCells(move)) {
			Cell c = move[ic];
			score -= 100 * c.tower * c.tower;
			for (int id : CellsAround(ic)) {
				Cell d = move[id];
				if (d.builder == ' ' && !d.dome) {
					if (d.tower == 3)
						score = -INF;
					score -= d.tower * d.tower;
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
		score += 100 * c.tower * c.tower;
		for (int id : CellsAround(ic)) {
			Cell d = move[id];
			if (d.builder == ' ' && !d.dome)
				score += d.tower * d.tower;
			// blocking with domes
			Cell d2 = state[id];
			if (d.dome && !d2.dome && d.tower == 3) {
				bool blocking = false;
				for (int ie : CellsAround(id)) {
					Cell e = move[ie];
					if (e.builder != ' ' && e.tower == 2 && Player(e) != state.player)
						blocking = true;
				}
				if (blocking)
					score += 100000;
			}
		}
	}
	for (int ic : PlayerCells(move)) {
		Cell c = move[ic];
		score -= 100 * c.tower * c.tower;
		for (int id : CellsAround(ic)) {
			Cell d = move[id];
			if (d.builder == ' ' && !d.dome) {
				if (d.tower == 3)
					score = -INF;
				score -= d.tower * d.tower;
			}
		}
	}
	return score;
}

double SubTreeScore(const State& prev, const State& state, int depth, bool maxi) {
	if (depth == 0)
		return Heuristic(prev, state) * (maxi ? 1 : -1);

	vector<State> moves;
	GenerateTurns(state, moves);
	if (moves.size() == 0)
		return maxi ? -INF : INF;
	if (moves.size() == 1 && moves[0].victory)
		return maxi ? INF : -INF;

	double best = (depth % 2) ? INF : -INF;
	for (const State& m : moves) {
		double s = SubTreeScore(state, m, depth - 1, !maxi);
		best = maxi ? max(best, s) : min(best, s);
	}
	return best;
}

optional<State> BruteBot(const State& state) {
	vector<State> moves;
	GenerateTurns(state, moves);
	if (moves.size() == 0)
		return nullopt;
	if (moves.size() == 1 && moves[0].victory)
		return moves[0];

	return bestElement(cspan<State>(moves), [&state](const State& m) {
		return SubTreeScore(state, m, 2, false);
	});
}

using Strategy = std::function<optional<State>(const State&)>;

double RelativeSkill2(int games, vector<God> gods, Strategy a, Strategy b) {
	vector<Strategy> strategy = { a, b };
	vector<int> wins;
	wins.resize(gods.size(), 0);
	for (int i = 0; i < games; i++) {
		State state;
		PlaceBuilder(state, 'A');
		PlaceBuilder(state, 'a');
		PlaceBuilder(state, 'B');
		PlaceBuilder(state, 'b');
		for (int j = 0; j < state.gods.size(); j++)
			state.gods[j] = (j < gods.size()) ? gods[j] : God::Dead;

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

double RelativeSkill(int games, vector<God> gods, Strategy a, Strategy b) {
	return (RelativeSkill2(games / 2, gods, a, b) + 1 - RelativeSkill2(games / 2, {gods[1], gods[0]}, b, a)) / 2;
}

void SingleMatch(vector<God> gods, Strategy a, Strategy b) {
	State state;
	PlaceBuilder(state, 'A');
	PlaceBuilder(state, 'a');
	PlaceBuilder(state, 'B');
	PlaceBuilder(state, 'b');
	for (int j = 0; j < state.gods.size(); j++)
		state.gods[j] = (j < gods.size()) ? gods[j] : God::Dead;
	Print(state);

	vector<Strategy> strategy = { a, b };
	while (!state.victory) {
		auto s = strategy[state.player](state);
		if (!s.has_value())
			break;
		state = *s;
		Print(state);
	}
}

int main(int argc, char* argv[]) {
	srand(time(0));
	vector<God> gods = { God::None, God::None };

	print("Greedy - Rand %s\n", RelativeSkill(10000, gods, GreedyBot, RandBot));
	print("Brute - Rand %s\n", RelativeSkill(100, gods, BruteBot, RandBot));
	print("Brute - Greedy %s\n", RelativeSkill(100, gods, BruteBot, GreedyBot));
	return 0;
}
