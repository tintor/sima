#ifndef __UTIL_H__
#define __UTIL_H__

#include <strstream>
#include <algorithm>

// TODO implement without a swap
template<typename T>
T median(T x, T y, T z) {
    if (x > y) std::swap(x, y);
    if (y > z) std::swap(y, z);
    if (x > y) std::swap(x, y);
    return y;
}

template<typename Container>
void release(Container& container) {
	Container empty;
	std::swap(empty, container);
}

template<typename Container>
void sort(Container& container) {
	std::sort(container.begin(), container.end());
}

template<typename Container, typename Func>
void sort(Container& container, const Func& func) {
	std::sort(container.begin(), container.end(), func);
}

template<typename T> constexpr
T min(T a, T b, T c) { return std::min(a, std::min(b, c)); }
template<typename T> constexpr
T max(T a, T b, T c) { return std::max(a, std::max(b, c)); }
template<typename T> constexpr
T min(T a, T b, T c, T d) { return std::min(std::min(a, b), std::min(c, d)); }
template<typename T> constexpr
T max(T a, T b, T c, T d) { return std::max(std::max(a, b), std::max(c, d)); }

// from C++17, but not in OSX clang yet
template<typename Number>
void from_chars(const char* begin, const char* end, Number& out) {
    std::istrstream is(begin, end - begin);
    is >> out;
}

template<typename Vector>
void remove_dups(Vector& v) {
	std::sort(v.begin(), v.end());
	v.erase(std::unique(v.begin(), v.end()), v.end());
}

#endif
