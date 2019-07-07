#pragma once
#include "core/std.h"
#include "core/array_bool.h"
#include "core/murmur3.h"

using Agent = uint;
using Boxes = array_bool<32 * 3>;

// sizeof 16 bytes
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

// sizeof 8 bytes
struct StateInfo {
	ushort distance = 0; // pushes so far
	ushort heuristic = 0; // estimated pushes remaining
	short dir = -1;
	short prev_agent = -1;
};
