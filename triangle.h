#pragma once
#include "segment.h"
#include "range.h"
#include "scalar.h"
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
	format_s(s, "(%s, %s, %s)", v.a, v.b, v.c);
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

inline std::string wkt(const ipolygon2& poly) {
	std::string s;
	s += "POLYGON (";
	if (poly.size() > 0) {
		for (auto p : poly) {
			format_e(s, "", p);
			s += ", ";
		}
		format_e(s, "", poly.front());
	}
	s += ')';
	return s;
}

inline std::string wkt(const imesh2& mesh) {
	std::string s;
	s += "MULTIPOLYGON ((";
	for (auto i : range(mesh.size())) {
		if (i != 0)
			s += ", ";
		const itriangle2& m = mesh[i];
		s += '(';
		format_e(s, "", m.a);
		s += ", ";
		format_e(s, "", m.b);
		s += ", ";
		format_e(s, "", m.c);
		s += ", ";
		format_e(s, "", m.a);
		s += ')';
	}
	s += "))";
	return s;
}

// overflow safe
inline long edge_area(ivec2 a, ivec2 b) {
	return ((long)a.x + (long)b.x) * ((long)a.y - (long)b.y);
}

// throws on overflow
inline long area(const ipolygon2& poly) {
	long area = 0;
	if (poly.size() > 0) {
	    auto a = poly.back();
		for (auto b : poly) {
			// This is a little faster than calling edge_area instead?
			area = addi(area, muli(addi(a.x, b.x), subi(a.y, b.y)));
			//area = addi(area, edge_area(a, b));
			a = b;
		}
	}
	return area;
}

// throws on overflow
inline long area(ivec2 a, ivec2 b, ivec2 c) {
	long ab = edge_area(a, b);
	long bc = edge_area(b, c);
	long ca = edge_area(c, a);
	long z = addi(ab, bc);
	return addi(z, ca);
}

// throws on overflow
inline long area(itriangle2 m) {
	return area(m.a, m.b, m.c);
}
