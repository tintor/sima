#pragma once

#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>

using std::thread;
using std::atomic;
using std::unique_lock;
using std::lock_guard;
using std::mutex;
using std::condition_variable;

inline void parallel(size_t max_threads, const std::function<void()>& func) {
	vector<thread> workers;
	auto m = max_threads;
	workers.reserve(m);
	for (decltype(m) i = 0; i < m; i++)
		workers.emplace_back(func);
	for (thread& w : workers)
		w.join();
}

inline void parallel(const std::function<void()>& func) {
	parallel(thread::hardware_concurrency(), func);
}

inline void parallel_for(size_t count, size_t max_threads, const std::function<void(size_t)>& func) {
	atomic<size_t> next = 0;
	parallel(max_threads, [count, &next, func]() {
		for (size_t task = next++; task < count; task = next++)
			func(task);
	});
}

inline void parallel_for(size_t count, const std::function<void(size_t)>& func) {
	parallel_for(count, thread::hardware_concurrency(), func);
}
