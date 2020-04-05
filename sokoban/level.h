#pragma once
#include "core/std.h"
#include "core/small_bfs.h"
#include "core/thread.h"
#include "core/format.h"
#include "core/exception.h"
#include "phmap/phmap.h"

#include "sokoban/code.h"
#include "sokoban/cell.h"
#include "sokoban/state.h"
#include "sokoban/util.h"

struct Level {
	string name;
	int width; // for printing only
	vector<char> buffer; // xy -> code, for printing only

	// TODO unique_ptr<Cell>
	vector<Cell*> cells; // ordinal -> Cell*, goals first, then alive, then dead cells
	cspan<Cell*> alive() const { return cspan<Cell*>(cells.data(), num_alive); }
	cspan<Cell*> goals() const { return cspan<Cell*>(cells.data(), num_goals); }

	int num_goals;
	int num_alive;
	int num_boxes;

	DynamicState start;
};

inline Cell* GetCell(const Level* level, uint xy) {
	for (Cell* c : level->cells)
		if (c->xy == xy)
			return c;
	THROW(runtime_error, "xy = %s", xy);
}

template<typename Boxes>
string_view Emoji(
		const Level* level,
		Agent agent,
		const Boxes& boxes,
		uint xy,
		const Boxes& frozen,
		std::function<string_view(Cell*)> fn) {
	if (level->buffer[xy] == Code::Wall)
		return "✴️ ";
	if (level->buffer[xy] == 'e')
		return "  ";

	Cell* c = GetCell(level, xy);
	string_view e = fn(c);
	if (e != "")
		return e;

	if (c->id == agent)
		return c->goal ? "😎" : "😀";
	if (!c->alive)
	 	return "🌀";
	if (boxes[c->id]) {
		if (c->goal)
			return frozen[c->id] ? "Ⓜ️ " : "🔵";
		return "🔴";
	}
	if (c->goal)
		return "🏳 ";
	return "🕸️ ";
}

template<typename State>
void Print(const Level* level, const State& key, std::function<string_view(Cell*)> fn = []LAMBDA("")) {
	small_bfs<const Cell*> visitor(level->cells.size());
	auto frozen = goals_with_frozen_boxes(level->cells[key.agent], key.boxes, level->goals(), visitor);
	for (uint xy = 0; xy < level->buffer.size(); xy++) {
		print("%s", Emoji(level, key.agent, key.boxes, xy, frozen, fn));
		if (xy % level->width == level->width - 1)
			print("\n");
	}
}

int NumberOfLevels(string_view filename);
const Level* LoadLevel(string_view filename);
uint CellCount(string_view filename);
void PrintInfo(const Level* level);
