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

#endif
