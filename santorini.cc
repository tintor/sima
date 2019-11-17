#include <core/format.h>
#include <core/each.h>
#include <core/dynamic_array.h>
#include <core/exception.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct Cell {
	bool dome = false;
	char tower = 0; // 0-3 height of tower
	char builder = ' '; // Aa - first player, Bb - second player, ... (uppercase male, lowercase female)
};

struct State {
	bool canMoveUp = true; // Athena's power
	int player = 0; // player about to play
	array<array<Cell, 5>, 5> cell;
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

auto PrevPlayerCells(const State& state) {
	static_vector<int, 4> data;
	for (int i = 0; i < 25; i++)
		if (Player(state[i]) == 1 - state.player) // TODO 2 player assumption
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

auto CellsAround(int a) {
	static_vector<int, 8> data;
	int row = a / 5;
	int col = a % 5;
	for (int dr = -1; dr <= 1; dr++)
		for (int dc = -1; dc <= 1; dc++)
			if (dr != 0 || dc != 0) {
				int r = row + dr;
				int c = col + dc;
				if (0 <= r && r < 5 && 0 <= c && c < 5)
					data.push_back(r * 5 + c);
			}
	return data;
}

enum class God { None, Apollo, Artemis, Athena, Atlas, Demeter, Hephaestus, Hermes, Minotaur, Pan, Prometheus };
// Gods implemented: Apollo, Athena, Atlas, Pan

void GenerateMoves(cspan<God> god, const State& state, vector<State>& moves) {
	const int p = state.player;
	for (int ic : PlayerCells(state)) {
		Cell c = state[ic];
		// perform all possible moves
		for (int id : CellsAround(ic)) {
			Cell d = state[id];
			bool no_builder = d.builder == ' ' || (Player(d) != p && god[p] == God::Apollo);
			if (no_builder && !d.dome && d.tower - c.tower <= (state.canMoveUp ? 1 : 0)) {
				State s = state;
				swap(s[ic].builder, s[id].builder);
				if (god[p] == God::Pan && d.tower - c.tower <= -2) {
					s.victory = true;
					moves.clear();
					moves.push_back(s);
					return;
				}
				if (d.tower == 3) {
					s.victory = true;
					moves.clear();
					moves.push_back(s);
					return;
				}

				s.player = (s.player + 1) % god.size();
				if (god[s.player] == God::Athena)
					s.canMoveUp = true;
				else if (god[p] == God::Athena)
					s.canMoveUp = d.tower - c.tower < 1;

				/*if (god[s.player] == God::Artemis) {
					// perform second optional move
					for (int ie : CellsAround(id)) {
						Cell e = s[ie];
						bool no_builder2 = d.builder == ' ';
						if (no_builder2 && !d.dome && d.tower - c.tower <= (state.canMoveUp ? 1 : 0)) {
							State s2 = s;
							swap(s2[ic].builder, s[id].builder);
							if (e.tower == 3) {
								s2.victory = true;
								moves.clear();
								moves.push_back(s2);
								return;
							}

							// TODO perform all possible builds
						}
					}
				}*/

				// perform all possible builds
				for (int ie : CellsAround(id)) {
					Cell e = s[ie];
					if (e.dome || e.builder != ' ')
						continue;

					if (e.tower == 3 || god[p] == God::Atlas) {
						State s2 = s;
						s2[ie].dome = true;
						moves.push_back(s2);
					}
					if (e.tower < 3) {
						State s2 = s;
						s2[ie].tower += 1;
						moves.push_back(s2);
					}
				}
			}
		}
	}
}

void Print(const State& state) {
	print("player %s", state.player);
	if (!state.canMoveUp)
		print(", can't move up");
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

optional<State> Human(cspan<God> players, const State& state) {
	vector<State> moves;
	GenerateMoves(players, state, moves);
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

optional<State> RandBot(cspan<God> players, const State& state) {
	vector<State> moves;
	GenerateMoves(players, state, moves);
	if (moves.size() == 0)
		return nullopt;
	return moves[rand() % moves.size()];
}

template<typename T, typename E>
T bestElement(cspan<T> data, const E& func) {
	vector<size_t> best = {0};
	auto bestScore = func(data[0]);
	for (size_t i = 1; i < data.size(); i++) {
		auto score = func(data[i]);
		if (score > bestScore) {
			bestScore = score;
			best = {i};
		} else if (score == bestScore)
			best.push_back(i);
	}
	return data[best[rand() % best.size()]];
}

optional<State> GreedyBot(cspan<God> players, const State& state) {
	vector<State> moves;
	GenerateMoves(players, state, moves);
	if (moves.size() == 0)
		return nullopt;

	return bestElement(cspan<State>(moves), [&state](const State& move){
		double score = 0;
		for (int ic : PrevPlayerCells(move)) {
			Cell c = move[ic];
			score += 100 * c.tower * c.tower;
			for (int id : CellsAround(ic)) {
				Cell d = state[id];
				if (d.builder == ' ' && !d.dome)
					score += d.tower * d.tower;
			}
		}
		return score;
	});
}

optional<State> GreedyBot2(cspan<God> players, const State& state) {
	vector<State> moves;
	GenerateMoves(players, state, moves);
	if (moves.size() == 0)
		return nullopt;

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
						if (e.builder != ' ' && Player(e) != state.player)
							blocking = true;
					}
					if (blocking)
						score += 10000;
				}
			}
		}
		return score;
	});
}

optional<State> BruteBot(cspan<God> players, const State& state) {
	// TODO
	return state;
}

using Strategy = std::function<optional<State>(cspan<God>, const State&)>;

double RelativeSkill2(int games, vector<God> players, Strategy a, Strategy b) {
	vector<Strategy> strategy = { a, b };
	vector<int> wins;
	wins.resize(players.size(), 0);
	for (int i = 0; i < games; i++) {
		State state;
		PlaceBuilder(state, 'A');
		PlaceBuilder(state, 'a');
		PlaceBuilder(state, 'B');
		PlaceBuilder(state, 'b');

		while (true) {
			auto s = strategy[state.player](players, state);
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

double RelativeSkill(int games, vector<God> players, Strategy a, Strategy b) {
	return (RelativeSkill2(games / 2, players, a, b) + 1 - RelativeSkill2(games / 2, {players[1], players[0]}, b, a)) / 2;
}

void SingleMatch(vector<God> players, Strategy a, Strategy b) {
	State state;
	PlaceBuilder(state, 'A');
	PlaceBuilder(state, 'a');
	PlaceBuilder(state, 'B');
	PlaceBuilder(state, 'b');
	Print(state);

	vector<Strategy> strategy = { a, b };
	while (!state.victory) {
		auto s = strategy[state.player](players, state);
		if (!s.has_value())
			break;
		state = *s;
		Print(state);
	}
}

int main(int argc, char* argv[]) {
	srand(time(0));
	vector<God> players = { God::Athena, God::Apollo };
	//SingleMatch(players, Human, GreedyBot);
	print("Greedy2 - Greedy2 %s\n", RelativeSkill(10000, players, GreedyBot2, GreedyBot2));

	/*print("Greedy - Rand %s\n", RelativeSkill(10000, players, GreedyBot, RandBot));
	print("Greedy2 - Rand %s\n", RelativeSkill(10000, players, GreedyBot2, RandBot));
	print("Greedy - Greedy2 %s\n", RelativeSkill(10000, players, GreedyBot, GreedyBot2));*/
	return 0;
}
