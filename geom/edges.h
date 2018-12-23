#pragma once
#include <core/std.h>
#include <geom/segment.h>
#include <core/each.h>
#include <core/auto.h>
#include <type_traits>

// Use to iterate over all segment<Vec> generated from triangle<Vec> or polygon<Vec>
template<typename Vec>
struct edges_flat {
	int vertex = 0;
	cspan<Vec> s;

	edges_flat(cspan<Vec> s) : s(s) { }
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
constexpr auto Edges(cspan<T> poly) {
	return iterable(edges_flat(poly));
}

template<typename T>
constexpr auto Edges(span<T> poly) {
	return iterable(edges_flat(poly));
}
