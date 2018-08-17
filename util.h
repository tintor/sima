#pragma once
#include "std.h"
#include <algorithm>

template<typename T>
inline int sign(T a) {
	return (0 < a) - (a < 0);
}

template<size_t S>
bool aligned(const void* ptr) {
	return (reinterpret_cast<size_t>(ptr) % S) == 0;
}

template<typename T>
inline bool ordered(T a, T b, T c) {
	return a <= b && b <= c;
}

template<typename T>
T min(T a, T b) {
	return (a < b) ? a : b;
}

template<typename T, typename... Args>
T min(T a, Args... args) {
	return min(a, min(args...));
}

template<typename T>
T max(T a, T b) {
	return (a > b) ? a : b;
}

template<typename T, typename... Args>
T max(T a, Args... args) {
	return max(a, max(args...));
}

template<typename T>
T median(T x, T y, T z) {
	if (x <= y) {
		if (z <= x)
			return x;
		if (y <= z)
			return y;
		return z;
	}
	// y < x
	if (z <= y)
		return y;
	if (x <= z)
		return x;
	return z;
}

template<typename Container>
void sort(Container& container) {
	std::sort(container.begin(), container.end());
}

template<typename Container, typename Func>
void sort(Container& container, const Func& func) {
	std::sort(container.begin(), container.end(), func);
}

template<typename Vector>
void remove_dups(Vector& v) {
       std::sort(v.begin(), v.end());
       v.erase(std::unique(v.begin(), v.end()), v.end());
}
