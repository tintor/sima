#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

using std::atomic;
using std::condition_variable;
using std::lock_guard;
using std::mutex;
using std::thread;
using std::unique_lock;

inline void parallel(size_t max_threads, const std::function<void()>& func) {
    vector<thread> workers;
    auto m = max_threads;
    workers.reserve(m);
    for (decltype(m) i = 0; i < m; i++) workers.emplace_back(func);
    for (thread& w : workers) w.join();
}

inline void parallel(size_t max_threads, const std::function<void(size_t)>& func) {
    vector<thread> workers;
    auto m = max_threads;
    workers.reserve(m);
    for (decltype(m) i = 0; i < m; i++) workers.emplace_back(func, i);
    for (thread& w : workers) w.join();
}

inline void parallel(const std::function<void()>& func) { parallel(thread::hardware_concurrency(), func); }

inline void parallel(const std::function<void(size_t)>& func) { parallel(thread::hardware_concurrency(), func); }

inline void parallel_for(size_t count, size_t max_threads, const std::function<void(size_t)>& func) {
    atomic<size_t> next = 0;
    parallel(max_threads, [count, &next, func]() {
        for (size_t task = next++; task < count; task = next++) func(task);
    });
}

// MapFn: Result()
// ReduceFn: Result(Result a, Result b)

template<typename Result, typename MapFn, typename ReduceFn>
inline Result parallel_map_reduce(const MapFn& map_fn, const ReduceFn& reduce_fn, bool linear_reduce) {
    std::mutex mutex;
    std::optional<Result> result;
    parallel([&]() {
        std::optional<Result> my_result = map_fn();

        std::unique_lock lock(mutex);
        if (linear_reduce) {
            result = result.has_value() ? reduce_fn(std::move(result.value()), std::move(my_result.value())) : my_result;
        } else {
            while (result.has_value()) {
                std::optional<Result> local_result;
                std::swap(local_result, result);

                // TODO make unlocking exception safe
                lock.unlock();
                my_result.value() = reduce_fn(std::move(my_result.value()), std::move(local_result.value()));
                lock.lock();
            }
            result = std::move(my_result);
        }
    });
    return result.value();
}

inline void parallel_for(size_t count, const std::function<void(size_t)>& func) {
    parallel_for(count, thread::hardware_concurrency(), func);
}
