#pragma once
#include "std.h"
#include "segment.h"
#include "each.h"
#include "auto.h"
#include <type_traits>

// Use to iterate over all segment<Vec> generated from triangle<Vec> or polygon<Vec>
template<typename Vec>
struct edges_flat {
	int vertex = 0;
	span<const Vec> s;

	edges_flat(span<const Vec> s) : s(s) { }
	optional<segment<Vec>> next() {
		if (vertex >= s.size()) return nullopt;
		ON_SCOPE_EXIT(vertex += 1);
		return (vertex == 0) ? segment(s.back(), s[0]) : segment(s[vertex - 1], s[vertex]);
	}
};

template<typename T>
constexpr auto Edges(const vector<T>& poly) {
	return iterable(edges_flat<T>(poly));
}

template<typename T>
constexpr auto Edges(span<const T> poly) {
	return iterable(edges_flat(poly));
}

template<typename T>
constexpr auto Edges(span<T> poly) {
	return iterable(edges_flat(poly));
}
