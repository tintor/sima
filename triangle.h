#pragma once
#include "segment.h"
#include "range.h"
#include "scalar.h"
#include "edges.h"
#include <vector>

template<typename Vec>
struct triangle {
	Vec a, b, c;

	triangle() { }
	constexpr triangle(Vec a, Vec b, Vec c) : a(a), b(b), c(c) { }
	constexpr triangle(const triangle& v) : a(v.a), b(v.b), c(v.c) { }

	Vec& operator[](int idx) { return (&a)[idx]; }
	Vec operator[](int idx) const { return (&a)[idx]; }

	array<segment<Vec>, 3> edges() const { return {segment{a, b}, segment{b, c}, segment{c, a}}; }
	array<segment<Vec&>, 3> edges() { return {segment{a, b}, segment{b, c}, segment{c, a}}; }

	bool operator==(triangle v) const { return equal(a, v.a) && equal(b, v.b) && equal(c, v.c); }
	bool operator!=(triangle v) const { return !operator==(v); }

	auto begin() { return &a; }
	auto end() { return &a + 3; }

	using const_iterator = const Vec*;
	auto begin() const { return &a; }
	auto end() const { return &a + 3; }
};

using triangle2 = triangle<double2>;
using triangle3 = triangle<double3>;

template<typename Vec>
void format_e(string& s, string_view spec, triangle<Vec> v) {
	format_s(s, "(%s, %s, %s)", v.a, v.b, v.c);
}

// TODO move polygon and mesh stuff outside

using polygon2 = vector<double2>;
using polygon3 = vector<double3>;

// TODO remove
template<typename T>
inline bool equal(const vector<T>& a, const vector<T>& b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); i++)
		if (!equal(a[i], b[i]))
			return false;
	return true;
}

inline bool operator==(const polygon2& a, const polygon2& b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); i++)
		if (!equal(a[i], b[i]))
			return false;
	return true;
}

inline bool operator==(const polygon3& a, const polygon3& b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); i++)
		if (!equal(a[i], b[i]))
			return false;
	return true;
}

inline bool operator!=(const polygon3& a, const polygon3& b) {
	return !operator==(a, b);
}

using mesh2 = vector<triangle2>;
using mesh3 = vector<triangle3>;

inline string wkt(const polygon2& poly) {
	string s;
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

inline string wkt(const mesh2& mesh) {
	string s;
	s += "MULTIPOLYGON ((";
	for (auto i : range(mesh.size())) {
		if (i != 0)
			s += ", ";
		const triangle2& m = mesh[i];
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

inline double edge_area(double2 a, double2 b) {
	return (a.x + b.x) * (a.y - b.y);
}

inline double area(const polygon2& poly) {
	double area = 0;
	if (poly.size() > 0)
		for (auto [a, b] : edgesOf(poly))
			area += edge_area(a, b);
	return area;
}

inline double area(double2 a, double2 b, double2 c) {
	return edge_area(a, b) + edge_area(b, c) + edge_area(c, a);
}

inline double area(triangle2 m) {
	return area(m.a, m.b, m.c);
}

inline double3 compute_normal(triangle3 v) { return compute_normal(v.a, v.b, v.c); }
