#pragma once
#include "core/std.h"

#include "sokoban/level.h"
#include "sokoban/frozen.h"
#include "sokoban/util.h"

using Corral = vector<uchar>; // avoid slower bit optimized vector<bool>

bool is_single_component(const Level* level, const Corral& corral);

inline void add(vector<uchar>& dest, const vector<uchar>& src) {
	for (size_t i = 0; i < min(dest.size(), src.size()); i++)
		if (src[i])
			dest[i] = true;
}

template<typename State>
class Corrals {
public:
	using Boxes = typename State::Boxes;

	Corrals(const Level* level) : _level(level), _corral(level->cells.size()), _reachable(level->cells.size()) {}
	void find_unsolved_picorral(const State& s);

	bool has_picorral() const { return _has_picorral; }
	const Corral& picorral() const { return _picorral; }
	optional<Corral> opt_picorral() const { if (_has_picorral) return _picorral; return nullopt; }

private:
	void find_corrals(const State& s);
	void add_if_picorral(const Boxes& boxes);

	const Level* _level;
	Corral _corral;
	vector<pair<Corral, bool>> _corrals;
	int _corrals_size;
	vector<uchar> _reachable;

	bool _has_picorral;
	int _picorral_pushes;
	Corral _picorral;
};

// corral with a goal without a box OR a box not on goal
template<typename Boxes>
bool is_unsolved_corral(const Level* level, const Boxes& boxes, const Corral& corral) {
	for (Cell* a : level->alive())
		if (corral[a->id] && a->goal != boxes[a->id])
			return true;
	return false;
}

// I-Corral is a corral where all boxes on the barrier can only be pushed inwards (into the corral) by the first push
// PI-Corral is an I-Corral where the player can perform all legal first pushes into the corral,
// meaning the player can reach all the relevant boxes from all relevant directions.
template<typename Boxes>
bool is_picorral(const Level* level, const Boxes& boxes, const vector<uchar>& reachable, const Corral& corral, int& count) {
	for (Cell* a : level->alive())
		if (corral[a->id] && boxes[a->id])
			for (auto [b, q] : a->pushes) {
				if (!corral[b->id] && !corral[q->id])
					return false;
				if (!boxes[b->id] && corral[b->id] && !corral[q->id]) {
					count += 1;
					if (boxes[q->id]) {
						if (is_frozen_on_goal_simple(q, boxes))
							continue;
						return false;
					}
					Boxes nboxes(boxes);
					nboxes.reset(a->id);
					nboxes.set(b->id);
					if (is_simple_deadlock(b, nboxes))
						continue;
					if (!reachable[q->id])
						return false;
				}
			}
	return true;
}

template<typename State>
void Corrals<State>::find_corrals(const State& s) {
	// optimize: memory allocation
	small_bfs<Cell*> visitor(_level->cells.size());
	visitor.clear();
	visitor.add(_level->cells[s.agent], s.agent);
	for (auto& e : _reachable)
		e = false;
	for (Cell* a : visitor) {
		_reachable[a->id] = true;
		for (auto [_, b] : a->moves)
			if (!s.boxes[b->id])
				visitor.add(b, b->id);
			else
				_reachable[b->id] = true;
	}

	_corrals.clear();
	for (Cell* q : _level->cells)
		if (!s.boxes[q->id] && !visitor.visited[q->id]) {
			// optimize: memory allocation
			vector<uchar> corral(_level->cells.size(), false);

			visitor.add(q, q->id);
			for (Cell* a : visitor) {
				corral[a->id] = true;

				for (Cell* b : a->dir8)
					if (b && !corral[b->id] && s.boxes[b->id])
						corral[b->id] = true;

				for (auto [_, b] : a->moves)
					if (!s.boxes[b->id])
						visitor.add(b, b->id);
					else
					corral[b->id] = true;
			}
			bool unsolved = is_unsolved_corral(_level, s.boxes, corral);
			_corrals.emplace_back(std::move(corral), unsolved);
		}
}

template<typename State>
void Corrals<State>::add_if_picorral(const Boxes& boxes) {
	int pushes = 0;
	if (is_picorral(_level, boxes, _reachable, _corral, /*out*/pushes))
		if (!_has_picorral || pushes < _picorral_pushes) {
			_picorral = _corral;
			_picorral_pushes = pushes;
			_has_picorral = true;
		}
}

// optimize: find the most expensive cases and improve them
template<typename State>
void Corrals<State>::find_unsolved_picorral(const State& s) {
	find_corrals(s);

	_has_picorral = false;
	_picorral_pushes = std::numeric_limits<int>::max();
	if (_corrals.size() >= 8) {
		// size 1
		for (const auto& c : _corrals)
			if (c.second) {
				_corral = c.first;
				add_if_picorral(s.boxes);
			}
		// size 2
		for (size_t a = 0; a < _corrals.size(); a++)
			for (size_t b = a + 1; b < _corrals.size(); b++)
				if (_corrals[a].second || _corrals[b].second) {
					// generate subset from individual corrals
					_corral = _corrals[a].first;
					add(_corral, _corrals[b].first);
					add_if_picorral(s.boxes);
				}
		// size all
		for (auto& e : _corral) e = false;
		for (size_t i = 0; i < _corrals.size(); i++)
			add(_corral, _corrals[i].first);
		if (is_unsolved_corral(_level, s.boxes, _corral))
			add_if_picorral(s.boxes);
	} else {
		for (size_t subset = 1; subset < (1lu << _corrals.size()); subset++) {
			// subset must contain at least one unsolved corral
			bool unsolved = false;
			for (size_t i = 0; i < _corrals.size(); i++) {
				size_t m = 1lu << i;
				if ((subset & m) == m && _corrals[i].second) {
					unsolved = true;
					break;
				}
			}
			if (!unsolved)
				continue;

			for (auto& e : _corral) e = false;
			for (size_t i = 0; i < _corrals.size(); i++) {
				size_t m = 1lu << i;
				if ((subset & m) == m)
					add(_corral, _corrals[i].first);
			}
			add_if_picorral(s.boxes);
		}
	}
}
