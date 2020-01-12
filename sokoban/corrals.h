#pragma once
#include "core/dynamic_array.h"

#include "sokoban/level.h"
#include "sokoban/frozen.h"
#include "sokoban/util.h"

using Corral = dynamic_array<bool>;

class Corrals {
public:
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
	dynamic_array<bool> _reachable;

	bool _has_picorral;
	int _picorral_pushes;
	Corral _picorral;
};
