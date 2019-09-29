#include "core/util.h"
#include "core/bits_util.h"
#include "core/format.h"
#include "core/small_bfs.h"
#include "core/array_deque.h"
#include "core/timestamp.h"
#include "core/thread.h"

#include "sokoban/level.h"
#include "sokoban/util.h"
#include "sokoban/hungarian.h"

#include "phmap/phmap.h"
using namespace std::chrono_literals;

struct StateMap {
	constexpr static int SHARDS = 32;
	array<mutex, SHARDS> locks;
	array<phmap::flat_hash_map<State, StateInfo>, SHARDS> data;

	static int shard(const State& s) {
		return fmix64(std::hash<Boxes>()(s.boxes) * 7) % SHARDS;
	}

	void lock(int shard) {
		Timestamp lock_ts;
		locks[shard].lock();
		overhead += lock_ts.elapsed();
	}

	void lock2(int shard) {
		Timestamp lock_ts;
		locks[shard].lock();
		overhead2 += lock_ts.elapsed();
	}

	void unlock(int shard) {
		locks[shard].unlock();
	}

	bool contains(const State& s, int shard) {
		return data[shard].find(s) != data[shard].end();
	}

	StateInfo get(const State& s, int shard) {
		return data[shard][s];
	}

	StateInfo* query(const State& s, int shard) {
		auto& d = data[shard];
		auto it = d.find(s);
		if (it == d.end())
			return nullptr;
		return &it->second;
	}

	void add(const State& s, const StateInfo& si, int shard) {
		data[shard].emplace(s, si);
	}

	long size() {
		long result = 0;
		for (int i = 0; i < SHARDS; i++) {
			locks[i].lock();
			result += data[i].size();
			locks[i].unlock();
		}
		return result;
	}

	atomic<long> overhead = 0;
	atomic<long> overhead2 = 0;
};

template<typename T>
void ensure_size(vector<T>& vec, size_t s) {
	if (s > vec.size())
		vec.resize(round_up_power2(s));
}

class StateQueue {
public:
	StateQueue() {
		queue.resize(256);
	}

	void push(const State& s, uint priority) {
		Timestamp lock_ts;
		queue_lock.lock();
		_overhead += lock_ts.elapsed();

		ensure_size(queue, priority + 1);
		queue[priority].push_back(s);
		if (priority < min_queue)
			min_queue = priority;
		bool notify = queue_size == 0;
		queue_size += 1;
		queue_lock.unlock();

		if (notify)
			queue_push_cv.notify_one();
	}

	optional<State> pop() {
		Timestamp lock_ts;
		std::unique_lock<mutex> lk(queue_lock);
		_overhead2 += lock_ts.elapsed();

		if (queue_size == 0) {
			blocked_on_queue += 1;
			while (queue_size == 0) {
				if (!running)
					return nullopt;
				if (blocked_on_queue >= thread::hardware_concurrency()) { // TODO remove this dependency
					running = false;
					queue_push_cv.notify_all();
					return nullopt;
				}
				queue_push_cv.wait(lk);
			}
			blocked_on_queue -= 1;
		}

		if (!running)
			return nullopt;
		while (queue[min_queue].size() == 0)
			min_queue += 1;
		State s = queue[min_queue].pop_front();
		queue_size -= 1;
		return s;
	}

	void shutdown() {
		std::unique_lock<mutex> lk(queue_lock);
		running = false;
		queue_push_cv.notify_all();
	}

	template<class Rep, class Period>
	bool wait_while_running_for(const std::chrono::duration<Rep, Period>& rel_time) {
		std::unique_lock<mutex> lk(queue_lock);
		if (running)
			queue_push_cv.wait_for(lk, rel_time);
		return running;
	}

	double overhead() {
		std::unique_lock<mutex> lk(queue_lock);
		return Timestamp::to_s(_overhead);
	}

	double overhead2() {
		std::unique_lock<mutex> lk(queue_lock);
		return Timestamp::to_s(_overhead2);
	}

private:
	bool running = true;
	long _overhead = 0;
	long _overhead2 = 0;
	uint min_queue = 0;
	uint blocked_on_queue = 0;
	uint queue_size = 0;
	vector<array_deque<State>> queue;
	mutex queue_lock;
	condition_variable queue_push_cv;
};

struct Solver {
	const Level* level;
	StateMap states;

	StateQueue queue;
	atomic<size_t> queue_pop_misses = 0;

	Solver(const Level* level)
		: level(level)
	{ }

	optional<pair<State, StateInfo>> queue_pop() {
		while (true) {
			optional<State> s = queue.pop();
			if (!s)
				return nullopt;

			int shard = StateMap::shard(*s);
			states.lock2(shard);
			StateInfo* q = states.query(*s, shard);
			if (!q && q->closed)  {
				states.unlock(shard);
				queue_pop_misses += 1;
				continue;
			}
			StateInfo si = *q; // copy before unlock
			q->closed = true;
			states.unlock(shard);
			return pair<State, StateInfo>{ *s, si };
		}
	}

	void normalize(State& s, small_bfs<Cell*>& visitor) {
		visitor.clear();
		visitor.add(level->cells[s.agent], s.agent);
		for (Cell* a : visitor)
			for (auto [_, b] : a->moves)
				if (!s.boxes[b->id]) {
					visitor.add(b, b->id);
					if (b->id < s.agent)
						s.agent = b->id;
				}
	}

	// +1 for every box not on goal
	uint heuristic(const State& s) {
		uint cost = 0;
		for (uint i = 0; i < level->num_goals; i++)
			if (!s.boxes[i])
				cost += 1;
		return cost;
	}

	optional<pair<State, StateInfo>> solve(int verbosity = 0) {
		State start = level->start;

		{
			small_bfs<Cell*> visitor(level->cells.size());
			normalize(start, visitor);
		}
		states.add(start, StateInfo(), StateMap::shard(start));
		queue.push(start, 0);

		Boxes goals;
		for (Cell* c : level->goals())
			goals[c->id] = true;

		if (start.boxes == goals)
			return pair<State, StateInfo>{ start, StateInfo() };

		optional<pair<State, StateInfo>> result;
		mutex result_lock;

		thread monitor([verbosity, this]() {
			Timestamp start_ts;
			while (verbosity > 0) {
				if (!queue.wait_while_running_for(5s))
					break;
				long seconds = std::lround(start_ts.elapsed_s());
				if (seconds < 4)
					continue;
				print("%s: states %dM, elapsed %d:%02d", level->name, states.size() / 1000000, seconds / 60, seconds % 60);
				print(", lock overhead (states1 %ds, states2 %ds", Timestamp::to_s(states.overhead.load()), Timestamp::to_s(states.overhead2.load()));
				print(", queue1 %ds, queue2 %ds)", queue.overhead(), queue.overhead2());
				print(", queue_pop_misses %s\n", queue_pop_misses.load());
				if (verbosity >= 2) {
					// TODO print level at the top of queue Print(level, s);
				}
			}
		});
		parallel([&]() {
			small_bfs<Cell*> agent_visitor(level->cells.size());
			small_bfs<Cell*> norm_visitor(level->cells.size());
			while (true) {
				auto p = queue_pop();
				if (!p)
					return;
				const State& s = p->first;
				const StateInfo& si = p->second;

				agent_visitor.clear();
				agent_visitor.add(level->cells[s.agent], s.agent);
				for (Cell* a : agent_visitor) {
					for (auto [d, b] : a->moves) {
						if (!s.boxes[b->id]) {
							agent_visitor.add(b, b->id);
							continue;
						}
						// push
						Cell* c = b->dir[d];
						if (!c || !c->alive || s.boxes[c->id])
							continue;
						State ns(b->id, s.boxes);
						ns.boxes[b->id] = false;
						ns.boxes[c->id] = true;
						if (is_simple_deadlock(c, ns.boxes))
							continue;
						//if (!is_reversible_push(ns, d, level) && contains_frozen_boxes(c, ns.boxes))
						//	continue;
						normalize(ns, norm_visitor);

						int shard = StateMap::shard(ns);
						states.lock(shard);

						StateInfo* q = states.query(ns, shard);
						if (q) {
							// existing state
							if (si.distance + 1 >= q->distance) {
								states.unlock(shard);
								continue;
							}
							q->dir = d;
							q->distance = si.distance + 1;
							// no need to update heuristic
							q->prev_agent = a->id;
							states.unlock(shard);
							queue.push(ns, uint(q->distance) + uint(q->heuristic));
							continue;
						}
						// new state
						StateInfo nsi;
						nsi.dir = d;
						nsi.distance = si.distance + 1;
						// TODO holding lock while computing heuristic!
						nsi.heuristic = heuristic(ns);
						nsi.prev_agent = a->id;

						if (ns.boxes == goals) {
							states.unlock(shard);
							queue.shutdown();
							lock_guard g(result_lock);
							if (!result)
								result = pair<State, StateInfo>{ ns, nsi };
							return;
						}

						states.add(ns, nsi, shard);
						states.unlock(shard);
						queue.push(ns, int(nsi.distance) + int(nsi.heuristic));
					}
				}
			}
		});
		monitor.join();
		return result;
	}
};

inline string cat(string_view s, const string& a) {
	return format("%s%s", s, a);
}

// Big improvements:
// - hungarian heuristic and over-estimation
// - frozen boxes deadlock and reversible pushes
// - deadlock pattern database
// - pi-corrals
// - tunnel macros

// Engineering:
// - parallel solve will be much more effective with more extensive deadlock checking AND more expensive heuristic
// - after certain size switch to dense compressed hash-map
// - restart search, but reuse deadlock pattern database from before -> results in smaller state hash-table

// - simple trick for symmetrical levels: expand only non-symmetrical nodes from root AND mark original root as closed AND mark all other non-symmetrical modes as closed
// - try it by manually modifying level

// - merge closed and open set into one with new color field in StateInfo (avoids removing and re-inserting on queue.pop)
// - measure lock contention
// - separate queue lock contention for push and for pull
// - maybe shard queue lock?

// timing 5:15s for microban2 without 132,105,109,130,115 with 16 shards (8 threads) - one level at a time
// timing 4:46s for microban2 without 132,105,109,130,115 with 64 shards (8 threads) - one level at a time

// - large queue locking overhead on microban1:144!

int main(int argc, char** argv) {
	InitSegvHandler();
	Timestamp::init();
	Timestamp start_ts;
	constexpr string_view prefix = "sokoban/levels/";
	constexpr string_view file = "original";

	atomic<int> total = 0;
	atomic<int> completed = 0;
	atomic<int> skipped = 0;
	auto num = NumberOfLevels(format("%s%s", prefix, file));
	vector<const Level*> levels;
	mutex levels_lock;

	parallel_for(num, 1, [&](size_t task) {
		string name = format("%s:%s%s%d", file, (task + 1) < 10 ? "0" : "", (task + 1) < 100 ? "0" : "", task + 1);
		if (name == "original:24") // syntax?
			return;
		if (name == "microban2:132") // large maze with one block
			return;
		if (name == "microban2:105" || name == "microban2:109") // >1 hour
			return;
		if (name == "microban2:130" || name == "microban2:115") // slow solves
			return;

		// microban 130:
		// baseline no locks:        closed 33865k, open 21668k, elapsed 5m0s
		// baseline 1 locked thread: closed 32694k, open 20951k, elapsed 5m0s
		//          8 locked thread: closed 42538k, open 28242k, elapsed 5m0s

		// microban 115:
		// 1 locked thread: closed 19019k, open 76819k, elapsed 5m0s
		// 8 locked thread: closed 23757k, open 94741k, elapsed 5m0s

		if (name == "microban1:044" || name == "microban1:144") {
			skipped += 1;
			return;
		}
		if (CellCount(cat(prefix, name)) > 256) {
			skipped += 1;
			return;
		}

		print("%s\n", name);
		auto level = LoadLevel(cat(prefix, name));
		levels_lock.lock();
		levels.push_back(level);
		levels_lock.unlock();
	});
	sort(levels, [](const Level* a, const Level* b) { return a->name < b->name; });

	parallel_for(levels.size(), 1, [&](size_t task) {
		auto level = levels[task];

		PrintInfo(level);
		total += 1;

		Solver solver(level);
		auto solved = solver.solve(2);
		if (solved) {
			completed += 1;
			print("%s: solved in %d pushes!\n", level->name, solved->second.distance);
			if (false) {
				small_bfs<Cell*> visitor(level->cells.size());
				State s = solved->first;
				StateInfo si = solved->second;
				while (si.distance > 0) {
					Print(level, s);
					State ps;
					ps.agent = si.prev_agent;
					ps.boxes = s.boxes;
					Cell* a = solver.level->cells[ps.agent];
					if (!a)
						THROW(runtime_error);
					Cell* b = a->dir[si.dir];
					if (!b || ps.boxes[b->id])
						THROW(runtime_error);
					Cell* c = b->dir[si.dir];
					if (!c || !ps.boxes[c->id])
						THROW(runtime_error);
					ps.boxes[c->id] = false;
					ps.boxes[b->id] = true;
					solver.normalize(ps, visitor);
					s = ps;
					si = solver.states.get(s, StateMap::shard(s));
				}
			}
		}
		else
			print("%s: no solution!\n", level->name);
	});

	int seconds = std::lround(start_ts.elapsed_s());
	print("solved %d/%d in %dm:%ds (skipped %d)\n", completed.load(), total.load(), seconds / 60, seconds % 60, skipped.load());
	return 0;
}
