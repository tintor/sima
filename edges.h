#pragma once
#include "std.h"
#include "segment.h"

// Use to iterate over all segment<Vec> generated from triangle<Vec> or polygon<Vec>
// TODO iterate over segment<Vec&>
template<typename Container>
class edgesOf {
public:
	constexpr edgesOf(const Container& cont) : _begin(cont.begin()), _end(cont.end()) {}

	struct iterator {
		typename Container::const_iterator _a, _b;
		auto operator*() { return segment{*_a, *_b}; }
		iterator& operator++() { ++_a; ++_b; return *this; }
		constexpr bool operator!=(iterator other) { return _b != other._b; }
	};

	constexpr iterator begin() const { auto e = _end; return iterator{--e, _begin}; }
	constexpr iterator end() const { return iterator{_end, _end}; }

private:
	typename Container::const_iterator _begin, _end;
};
