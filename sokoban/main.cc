#include "core/auto.h"
#include "core/util.h"
#include "core/string_util.h"
#include "core/bits_util.h"
#include "core/format.h"
#include "core/small_bfs.h"
#include "core/array_deque.h"
#include "core/timestamp.h"
#include "core/thread.h"

#include "sokoban/level.h"
#include "sokoban/frozen.h"
#include "sokoban/util.h"
#include "sokoban/hungarian.h"

#include "phmap/phmap.h"
using namespace std::chrono_literals;

struct StateMap {
	constexpr static int SHARDS = 64;
	array<mutex, SHARDS> locks;
	array<phmap::flat_hash_map<State, StateInfo>, SHARDS> data;

	static int shard(const State& s) {
		return fmix64(s.boxes.hash() * 7) % SHARDS;
	}

	void print_sizes() {
		print("states map");
		for (int i = 0; i < SHARDS; i++) {
			locks[i].lock();
			print(" %h", data[i].size());
			locks[i].unlock();
		}
		print("\n");
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
			queue_push_cv.notify_all();
	}

	optional<pair<State, uint>> top() {
		uint mq;
		auto s = pop(&mq);
		if (s.has_value())
			return pair<State, uint>{*s, mq};
		return nullopt;
	}

	optional<State> pop(uint* top_ptr = nullptr) {
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

		if (top_ptr) {
			*top_ptr = min_queue;
			return queue[min_queue][0];
		}
		State s = queue[min_queue].pop_front();
		queue_size -= 1;
		return s;
	}

	size_t size() {
		std::unique_lock<mutex> lk(queue_lock);
		return queue_size;
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

// check if unreachable area can only be reached by pushing a box into it
bool is_corral(const Level* level, const Boxes& boxes, const dynamic_array<bool>& reachable, const dynamic_array<bool>& corral) {
	bool r = false;
	for (uint i = 0; i < level->num_alive; i++)
		if (boxes[i] && corral[i]) {
			Cell* a = level->cells[i];
			for (auto [d, p] : a->moves) {
				Cell* q = a->dir(d ^ 2);
				if (q && p->alive && !boxes[q->id] && !boxes[p->id]) {
					if (!corral[q->id] && !corral[p->id])
						return false;
					if (reachable[q->id] && corral[p->id])
						r = true;
				}
			}
		}
	return r;
}

struct Counters {
	typedef ulong T;

	T simple_deadlocks = 0;
	T frozen_box_deadlocks = 0;
	T heuristic_deadlocks = 0;
	T corral_cuts = 0;
	T duplicates = 0;
	T updates = 0;

	T total_ticks = 0;
	T queue_pop_ticks = 0;
	T corral_ticks = 0;
	T corral2_ticks = 0;
	T is_simple_deadlock_ticks = 0;
	T is_reversible_push_ticks = 0;
	T contains_frozen_boxes_ticks = 0;
	T norm_ticks = 0;
	T states_query_ticks = 0;
	T heuristic_ticks = 0;
	T state_insert_ticks = 0;
	T queue_push_ticks = 0;

	Counters() {
		memset(this, 0, sizeof(Counters));
	}

	void add(const Counters& src) {
		const T* s = reinterpret_cast<const T*>(&src);
		const T* e = s + sizeof(Counters) / sizeof(T);
		static_assert(sizeof(Counters) == sizeof(T) * 18);
		T* d = reinterpret_cast<T*>(this);
		while (s != e)
			*d++ += *s++;
	}
};

struct Solver {
	const Level* level;
	StateMap states;
	StateQueue queue;
	vector<Counters> counters;

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
			if (!q || q->closed)  {
				states.unlock(shard);
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
		for (Cell* a : visitor) {
			if (a->id < s.agent)
				s.agent = a->id;
			for (auto [_, b] : a->moves)
				if (!s.boxes[b->id])
					visitor.add(b, b->id);
		}
	}

	uint heuristic_simple(const Boxes& boxes) {
		uint cost = 0;
		for (uint i = 0; i < level->num_alive; i++)
			if (boxes[i])
				cost += level->cells[i]->min_push_distance;
		return cost;
	}

	// excludes frozen goals from costs
	uint heuristic(const Boxes& boxes) {
		array<uint, Boxes::size()> goal;
		int goals = 0;
		for (Cell* g : level->goals())
			if (!boxes[g->id] || !is_frozen_on_goal_simple(g, boxes))
				goal[goals++] = g->id;
		if (goals == level->num_goals)
			return heuristic_simple(boxes);

		uint cost = 0;
		for (Cell* box : level->alive())
			if (boxes[box->id] && !box->goal) {
				// min push distance out of all non-frozen goals
				uint dist = Cell::Inf;
				for (uint i = 0; i < goals; i++)
					minimize(dist, box->push_distance[goal[i]]);
				cost += dist;
			}
		return cost;
	}

	optional<pair<State, StateInfo>> solve(int verbosity = 0) {
		Timestamp start_ts;
		State start = level->start;

		{
			small_bfs<Cell*> visitor(level->cells.size());
			normalize(start, visitor);
		}
		states.add(start, StateInfo(), StateMap::shard(start));
		queue.push(start, 0);

		Boxes goals;
		for (Cell* c : level->goals())
			goals.set(c->id);

		if (start.boxes == goals)
			return pair<State, StateInfo>{ start, StateInfo() };

		optional<pair<State, StateInfo>> result;
		mutex result_lock;

		counters.resize(thread::hardware_concurrency());
		thread monitor([verbosity, this, &result_lock, start_ts]() {
			dynamic_array<bool> reachable;
			reachable.resize(level->cells.size());
			dynamic_array<bool> corral;
			corral.resize(level->cells.size());

			small_bfs<Cell*> agent_visitor(level->cells.size());

			bool running = verbosity > 0;
			while (running) {
				if (!queue.wait_while_running_for(5s))
					running = false;
				long seconds = std::lround(start_ts.elapsed_s());
				if (seconds < 4 && running)
					continue;

				auto total = states.size();
				auto open = queue.size();
				auto closed = (total >= open) ? total - open : 0;

				lock_guard g(result_lock);
				print("%s: states %h (%h %h %.1f)\n", level->name, total, closed, open, 100. * open / total);

				Counters q;
				for (const Counters& c : counters)
					q.add(c);
				ulong else_ticks = q.total_ticks;
				else_ticks -= q.queue_pop_ticks;
				else_ticks -= q.corral_ticks;
				else_ticks -= q.corral2_ticks;
				else_ticks -= q.is_simple_deadlock_ticks;
				else_ticks -= q.is_reversible_push_ticks;
				else_ticks -= q.contains_frozen_boxes_ticks;
				else_ticks -= q.norm_ticks;
				else_ticks -= q.states_query_ticks;
				else_ticks -= q.heuristic_ticks;
				else_ticks -= q.state_insert_ticks;
				else_ticks -= q.queue_push_ticks;

				print("elapsed %t (%Y | queue_pop %Y, corral %Y %Y, "
					"is_simple_deadlock %Y, is_reversible_push %Y, contains_frozen_boxes %Y, "
					"norm %Y, states_query %Y, heuristic %Y, "
					"state_insert %Y, queue_push %Y, else %Y)\n",
					seconds, q.total_ticks, q.queue_pop_ticks, q.corral_ticks, q.corral2_ticks,
					q.is_simple_deadlock_ticks, q.is_reversible_push_ticks, q.contains_frozen_boxes_ticks,
					q.norm_ticks, q.states_query_ticks,
					q.heuristic_ticks, q.state_insert_ticks, q.queue_push_ticks, else_ticks);
				print("deadlocks (simple %h, frozen_box %h, heuristic %h)",
					q.simple_deadlocks, q.frozen_box_deadlocks, q.heuristic_deadlocks);
				print(", corral cuts %h, dups %h, updates %h", q.corral_cuts, q.duplicates, q.updates);
				print(", locks (%Y %Y", states.overhead, states.overhead2);
				print(" %t %t)\n", queue.overhead(), queue.overhead2());
				//states.print_sizes();
				if (verbosity >= 2) {
					auto ss = queue.top();
					if (ss.has_value()) {
						State s = ss->first;

						int shard = StateMap::shard(s);
						states.lock(shard);
						StateInfo* q = states.query(s, shard);
						states.unlock(shard);
						print("queue cost %s, distance %s, heuristic %s\n", ss->second, q->distance, q->heuristic);

						// TODO don't duplicate this logic
						// find all reachable cells
						agent_visitor.clear();
						agent_visitor.add(level->cells[s.agent], s.agent);
						reachable.fill(false);
						for (Cell* a : agent_visitor) {
							reachable[a->id] = true;
							for (auto [_, b] : a->moves)
								if (!s.boxes[b->id])
									agent_visitor.add(b, b->id);
						}

						// find all unreachable cells
						bool push_into_corral = false;
						for (Cell* a : level->cells)
							if (!s.boxes[a->id] && !agent_visitor.visited[a->id]) {
								corral.fill(false);
								bool candidate = false;

								agent_visitor.add(a, a->id);
								for (Cell* e : agent_visitor) {
									corral[e->id] = true;
									for (Cell* b : e->dir8) // TODO e->moves8
										if (b && s.boxes[b->id]) {
											corral[b->id] = true;
											if (!b->goal)
												candidate = true;
										}
									if (e->goal)
										candidate = true;

									for (auto [_, b] : e->moves)
										if (!s.boxes[b->id])
											agent_visitor.add(b, b->id);
								}

								if (candidate && is_corral(level, s.boxes, reachable, corral)) {
									push_into_corral = true;
									break;
								}
							}

						Print(level, s, [&](Cell* c) {
							if (!push_into_corral || !corral[c->id])
								return "";
							if (s.boxes[c->id])
								return c->goal ? "🔷" : "⚪";
							if (c->goal)
								return "❔";
							if (!c->alive)
								return "❕";
							return "▫️ ";
						});
					}
				}
			}
		});

		parallel([&](size_t thread_id) {
			Counters& q = counters[thread_id];
			small_bfs<Cell*> agent_visitor(level->cells.size());
			small_bfs<Cell*> tmp_visitor(level->cells.size());

			dynamic_array<bool> reachable;
			reachable.resize(level->cells.size());
			dynamic_array<bool> corral;
			corral.resize(level->cells.size());

			while (true) {
				Timestamp queue_pop_ts;
				ON_SCOPE_EXIT(q.total_ticks += queue_pop_ts.elapsed());

				auto p = queue_pop();
				if (!p)
					return;
				const State& s = p->first;
				const StateInfo& si = p->second;

				// TODO heuristic: order goals somehow (corner goals first) and try to find solution in goal order
				// -> if that doesn't work then do a regular search

				// corral = unreachable area surrounded by boxes which are either:
				// - frozen on goal OR (reachable by agent and only pushable inside)

				// if corral contains a goal (assuming num_boxes == num_goals) or one of its fence boxes isn't on goal
				// then that corral must be prioritized for push (ignoring all other corrals)

				Timestamp corral_ts;
				q.queue_pop_ticks += queue_pop_ts.elapsed(corral_ts);

				// find all reachable cells
				agent_visitor.clear();
				agent_visitor.add(level->cells[s.agent], s.agent);
				reachable.fill(false);
				for (Cell* a : agent_visitor) {
					reachable[a->id] = true;
					for (auto [_, b] : a->moves)
						if (!s.boxes[b->id])
							agent_visitor.add(b, b->id);
				}
				// TODO remember [a, d] pairs with boxes avoid generating it again later

				Timestamp corral2_ts;
				q.corral_ticks += corral_ts.elapsed(corral2_ts);

				// find all unreachable cells
				bool push_into_corral = false;
				for (Cell* a : level->cells)
					if (!s.boxes[a->id] && !agent_visitor.visited[a->id]) {
						corral.fill(false);
						bool candidate = false;

						agent_visitor.add(a, a->id);
						for (Cell* e : agent_visitor) {
							corral[e->id] = true;
							for (Cell* b : e->dir8) // TODO e->moves8
								if (b && s.boxes[b->id]) {
									corral[b->id] = true;
									if (!b->goal)
										candidate = true;
								}
							if (e->goal)
								candidate = true;

							for (auto [_, b] : e->moves)
								if (!s.boxes[b->id])
									agent_visitor.add(b, b->id);
						}

						if (candidate && is_corral(level, s.boxes, reachable, corral)) {
							push_into_corral = true;
							break;
						}
					}
				q.corral2_ticks += corral2_ts.elapsed();

				agent_visitor.clear();
				agent_visitor.add(level->cells[s.agent], s.agent);
				for (Cell* a : agent_visitor) {
					for (auto [d, b] : a->moves) {
						if (!s.boxes[b->id]) {
							agent_visitor.add(b, b->id);
							continue;
						}
						// push
						Cell* c = b->dir(d);
						if (!c || !c->alive || s.boxes[c->id])
							continue;
						if (push_into_corral && !corral[c->id]) {
							q.corral_cuts += 1;
							continue;
						}
						State ns(b->id, s.boxes);
						ns.boxes.reset(b->id);
						ns.boxes.set(c->id);

						Timestamp deadlock_ts;
						if (TIMER(is_simple_deadlock(c, ns.boxes), q.is_simple_deadlock_ticks)) {
							q.simple_deadlocks += 1;
							continue;
						}
						auto dd = d;
						auto bb = b;
						if (TIMER(!is_reversible_push(ns, dd, level, tmp_visitor), q.is_reversible_push_ticks)
								&& TIMER(contains_frozen_boxes(bb, ns.boxes, tmp_visitor), q.contains_frozen_boxes_ticks)) {
							q.frozen_box_deadlocks += 1;
							continue;
						}

						Timestamp norm_ts;
						normalize(ns, tmp_visitor);

						Timestamp states_query_ts;
						q.norm_ticks += norm_ts.elapsed(states_query_ts);

						int shard = StateMap::shard(ns);
						states.lock(shard);

						constexpr int Overestimate = 2;
						StateInfo* qs = states.query(ns, shard);
						if (qs) {
							q.duplicates += 1;
							if (si.distance + 1 >= qs->distance) {
								// existing state
								states.unlock(shard);
							} else {
								qs->dir = d;
								qs->distance = si.distance + 1;
								// no need to update heuristic
								qs->prev_agent = a->id;
								states.unlock(shard);
								queue.push(ns, uint(qs->distance) + uint(qs->heuristic) * Overestimate);
								q.updates += 1;
							}
							q.states_query_ticks += states_query_ts.elapsed();
							continue;
						}

						// new state
						StateInfo nsi;
						nsi.dir = d;
						nsi.distance = si.distance + 1;

						Timestamp heuristic_ts;
						q.states_query_ticks += states_query_ts.elapsed(heuristic_ts);

						// TODO holding lock while computing heuristic!
						uint h = heuristic(ns.boxes);
						q.heuristic_ticks += heuristic_ts.elapsed();

						if (h == Cell::Inf) {
							states.unlock(shard);
							q.heuristic_deadlocks += 1;
							continue;
						}
						nsi.heuristic = h;

						nsi.prev_agent = a->id;

						Timestamp state_insert_ts;
						states.add(ns, nsi, shard);
						states.unlock(shard);

						Timestamp queue_push_ts;
						q.state_insert_ticks += state_insert_ts.elapsed(queue_push_ts);

						queue.push(ns, int(nsi.distance) + int(nsi.heuristic) * Overestimate);
						q.queue_push_ticks += queue_push_ts.elapsed();

						if (ns.boxes == goals) {
							queue.shutdown();
							lock_guard g(result_lock);
							if (!result)
								result = pair<State, StateInfo>{ ns, nsi };
							return;
						}
					}
				}
			}
		});
		monitor.join();
		return result;
	}
};

const cspan<string_view> Blacklist = {
	"original:24", // syntax?
	"microban1:44",
	"microban2:132", // large maze with one block
	"microban3:47", // takes 2+ hours
	"microban4:75", // takes 1+ hours
	"microban4:85", // takes 1.5+ hours
	"microban4:92", // deadlocks easily
	"microban4:96", // takes .5+ hours
};

string solve(string_view file) {
	Timestamp start_ts;

	constexpr string_view prefix = "sokoban/levels/";

	atomic<int> total = 0;
	atomic<int> completed = 0;
	vector<const Level*> levels;
	mutex levels_lock;

	vector<string> skipped;
	vector<string> unsolved;

	if (file.find(":"sv) != string_view::npos) {
		auto level = LoadLevel(cat(prefix, file));
		if (level)
			levels.push_back(level);
		else
			skipped.emplace_back(split(file, {':', '/'}).back());
	} else {
		auto num = NumberOfLevels(cat(prefix, file));
		parallel_for(num, [&](size_t task) {
			string name = format("%s:%d", file, task + 1);
			if (contains(Blacklist, string_view(name)) || CellCount(cat(prefix, name)) > 256) {
				unique_lock g(levels_lock);
				skipped.emplace_back(split(name, {':', '/'}).back());
				return;
			}

			auto level = LoadLevel(cat(prefix, name));
			if (!level) {
				unique_lock g(levels_lock);
				skipped.emplace_back(split(name, {':', '/'}).back());
				return;
			}

			unique_lock g(levels_lock);
			levels.push_back(level);
		});
		sort(levels, [](const Level* a, const Level* b) { return !natural_less(a->name, b->name); });
	}

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
					Cell* b = a->dir(si.dir);
					if (!b || ps.boxes[b->id])
						THROW(runtime_error);
					Cell* c = b->dir(si.dir);
					if (!c || !ps.boxes[c->id])
						THROW(runtime_error);
					ps.boxes.reset(c->id);
					ps.boxes.set(b->id);
					solver.normalize(ps, visitor);
					s = ps;
					si = solver.states.get(s, StateMap::shard(s));
				}
			}
		} else {
			print("%s: no solution!\n", level->name);
			unique_lock g(levels_lock);
			unsolved.emplace_back(split(level->name, {':', '/'}).back());
		}
		print("\n");
	});

	sort(unsolved, natural_less);
	sort(skipped, natural_less);

	string result;
	format_s(result, "solved %d/%d in %T", completed, total, start_ts.elapsed());
	format_s(result, " unsolved %s", unsolved);
	format_s(result, " skipped %s", skipped);
	return result;
}

int main(int argc, char** argv) {
	InitSegvHandler();
	Timestamp::init();
	if (argc == 1) {
		vector<string> results;
		for (auto file : {"microban1", "microban2", "microban3", "microban4", "microban5"})
			results.emplace_back(solve(file));
		for (auto result : results)
			print("%s\n", result);
	}
	if (argc == 2)
		print("%s\n", solve(argv[1]));
	return 0;
}
