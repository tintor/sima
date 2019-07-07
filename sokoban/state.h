#pragma once
#include "core/std.h"
#include "core/array_bool.h"
#include "core/murmur3.h"

using Agent = uint;
using Boxes = array_bool<32 * 3>;

struct State {
	Boxes boxes;
	Agent agent;

	State() {}
	State(Agent agent, const Boxes& boxes) : agent(agent), boxes(boxes) {}
};

inline bool operator==(const State& a, const State& b) {
	return a.agent == b.agent && a.boxes == b.boxes;
}

namespace std {
	template<> struct hash<State> {
		size_t operator()(const State& a) const {
			return std::hash<Boxes>()(a.boxes) ^ fmix64(a.agent);
		}
	};
}

struct StateInfo {
	int dir = -1;
	uint pushes = 0;
	uint total_dist = 0;
	int prev_agent = -1;
	bool closed = false;
};
