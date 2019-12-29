#pragma once
#include <core/std.h>
#include <algorithm>

template<typename T>
inline int sign(T a) {
	return (0 < a) - (a < 0);
}

template<typename T>
T clamp(T a, T min, T max) {
	if (a > max) return max;
	if (a < min) return min;
	return a;
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
void minimize(T& a, T b) {
	if (b < a) {
		a = b;
	}
}

template<typename T>
void maximize(T& a, T b) {
	if (b > a) {
		a = b;
	}
}

template<typename T>
T min(T a, T b) {
	return (a < b) ? a : b;
}

template<typename T, typename... Args>
T min(T a, T b, Args... args) {
	return min(a, min(b, args...));
}

template<typename T>
T max(T a, T b) {
	return (a > b) ? a : b;
}

template<typename T, typename... Args>
T max(T a, T b, Args... args) {
	return max(a, max(b, args...));
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

template<typename Vector, typename Func>
void remove_dups(Vector& v, const Func& less, const Func& equal) {
       std::sort(v.begin(), v.end(), less);
       v.erase(std::unique(v.begin(), v.end(), equal), v.end());
}

template<typename T>
bool contains(const vector<T>& container, T element) {
	for (const T& e : container)
		if (e == element)
			return true;
	return false;
}

template<typename T>
bool contains(cspan<T> container, T element) {
	for (const T& e : container)
		if (e == element)
			return true;
	return false;
}

template<typename T, size_t M>
bool contains(const array<T, M>& container, T element) {
	return contains(cspan<T>(container), element);
}

template<typename Container, typename Func>
bool All(const Container& s, const Func& func) {
	for (const auto& e : s)
		if (!func(e))
			return false;
	return true;
}

template<typename Container, typename Func>
bool Any(const Container& s, const Func& func) {
	for (const auto& e : s)
		if (func(e))
			return true;
	return false;
}

template<typename T, typename P>
size_t IndexOfMax(cspan<T> s, P measure) {
	auto mv = measure(s[0]);
	int mi = 0;
	for (size_t i = 1; i < s.size(); i++) {
		auto v = measure(s[i]);
		if (v > mv) {
			mv = v;
			mi = i;
		}
	}
	return mi;
}
