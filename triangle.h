#pragma once
#include "segment.h"
#include <vector>

template<typename V>
struct triangle {
	V a, b, c;

	triangle() { }
	triangle(const triangle& v) : a(v.a), b(v.b), c(v.c) { }	
	triangle(V a, V b, V c) : a(a), b(b), c(c) { }	

	V& operator[](int idx) { return (&a)[idx]; }
	V operator[](int idx) const { return (&a)[idx]; }

	std::array<segment<V>, 3> edges() const { return {segment<V>{a, b}, segment<V>{b, c}, segment<V>{c, a}}; }
	std::array<segment<V&>, 3> edges() { return {segment<V&>{a, b}, {b, c}, {c, a}}; }

	bool operator==(const triangle& v) const { return a == v.a && b == v.b && c == v.c; }
	bool operator!=(const triangle& v) const { return a != v.a || b != v.b || c != v.c; }

	V* begin() { return &a; }
	V* end() { return &a + 3; }
	
	const V* begin() const { return &a; }
	const V* end() const { return &a + 3; }
};

template<typename T>
void format_e(std::string& s, std::string_view spec, triangle<T> v) {
	format_e(s, "", v.a);
	s += ", ";
	format_e(s, "", v.b);
	s += ", ";
	format_e(s, "", v.c);
}

using itriangle2 = triangle<ivec2>;
using itriangle3 = triangle<ivec3>;

using ipolygon2 = std::vector<ivec2>;
using ipolygon3 = std::vector<ivec3>;

// TODO move to own file
template<typename V>
struct polygon {
	std::vector<V> vertices;

	// TODO edges iterator
};

using imesh2 = std::vector<itriangle2>;
using imesh3 = std::vector<itriangle3>;
