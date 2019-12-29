#pragma once
#include "core/std.h"
#include "core/small_bfs.h"
#include "sokoban/state.h"

struct Level;

struct Cell {
	const Level* level;
	int id;
	int xy;
	bool goal;
	bool sink;
	bool alive;

	array<Cell*, 4> dir;
	dynamic_array<pair<int, Cell*>> moves;

	array<Cell*, 8> dir8;

	bool straight() const { return moves.size() == 2 && (moves[0].first ^ 2) == moves[1].first; }
};

struct Level {
	string name;
	int width; // for printing only
	vector<char> buffer; // xy -> code, for printing only

	vector<Cell*> cells; // ordinal -> Cell*, goals first, then alive, then dead cells
	auto alive() const { return span(cells.data(), num_alive); }
	auto goals() const { return span(cells.data(), num_goals); }

	int num_goals;
	int num_alive;
	int num_boxes;
	State start;
};

int NumberOfLevels(string_view filename);
const Level* LoadLevel(string_view filename);
uint CellCount(string_view filename);
void PrintInfo(const Level* level);
void Print(const Level* level, const State& key, std::function<string_view(Cell*)> fn = [](Cell*) {return "";});
