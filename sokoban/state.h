#pragma once
#include "core/std.h"
#include "core/array_bool.h"
#include "core/murmur3.h"

using Agent = uint;
using Boxes = array_bool<32 * 3>;

// TODO save space by combining State and StateInfo into one struct

// TODO save space by allocating minimal space for agent and boxes

// TODO save space by having separate maps for open and closed states

// TODO save space by keeping agent outside of state

// with State of 16 bytes and StateInfo of 8 bytes, solver lasted until 1265M states

// sizeof 16 bytes
struct State {
	Boxes boxes;
	Agent agent;

	State() {}
	State(Agent agent, const Boxes& boxes) : agent(agent), boxes(boxes) {}
};
static_assert(sizeof(State) == 16);

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

// TODO save space by not storing heuristic
// TODO save space by not storing distance (can be recomputed by backtracking)
// TODO remaining are:
//      - 1 bit - closed
//      - 2 bit - dir
//      - 7-8 bit - prev_agent

// original2 level with 70 cells, and 46 alive
// agent 7bit
// boxes 46bit
// closed 1bit
// dir 2bit
// prev_agent 7bit
// TOTAL 63 -> 8bytes!

// sizeof 8 bytes
struct StateInfo {
	ushort distance = 0; // pushes so far
	ushort heuristic = 0; // estimated pushes remaining
	char dir = -1;
	bool closed = false;
	short prev_agent = -1;
};
static_assert(sizeof(StateInfo) == 8);
