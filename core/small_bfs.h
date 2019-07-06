#pragma once
#include "core/each.h"
#include "core/exception.h"
#include "core/dynamic_array.h"

// very simple and fast array queue (but number of pushes is limited)
template<typename T>
class small_queue {
public:
	small_queue(uint capacity) : _data(capacity) { }
	bool empty() const { return _head == _tail; }
	operator bool() const { return !empty(); }
	void push(T e) { _data[_tail++] = e; }
	T pop() { return _data[_head++]; }
	void clear() { _head = _tail = 0; }
	uint tail() const { return _tail; }
private:
	dynamic_array<T> _data;
	uint _head = 0;
	uint _tail = 0;
};

template<typename T>
struct small_bfs : public each<small_bfs<T>> {
	small_queue<T> queue;
	dynamic_array<bool> visited;

	small_bfs(uint capacity) : queue(capacity) {
		visited.resize(capacity);
		memset(visited.begin(), 0, visited.size());
	}

	void clear() {
		queue.clear();
		memset(visited.begin(), 0, visited.size());
	}

	void add(T e, uint index) {
		if (!visited[index]) {
			visited[index] = true;
			queue.push(e);
		}
	}

	optional<T> next() {
		if (queue)
			return queue.pop();
		return nullopt;
	}
};
