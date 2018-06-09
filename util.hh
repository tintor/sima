#ifndef __UTIL_H__
#define __UTIL_H__

template<typename Container>
void release(Container& container) {
	Container empty;
	std::swap(empty, container);
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

#endif