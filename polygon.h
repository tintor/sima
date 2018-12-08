#pragma once
#include "triangle.h"

// TODO move polygon and mesh stuff outside

using polygon2 = vector<double2>;
using polygon3 = aligned_vector<double4>;

// polygon2 with holes
template<typename Vec>
class xpolygon {
public:
	xpolygon() {
		_offsets.push_back(0);
	}

	void add(Vec p) {
		assert(size() > 0);
		_vertices.push_back(p);
		_offsets.back() += 1;
	}

	void add(span<const Vec> poly) {
		_vertices << poly;
		_offsets.push_back(_offsets.back() + poly.size());
	}

	span<const Vec> operator[](uint idx) const {
		return span<const Vec>(_vertices.data() + _offsets[idx], _vertices.data() + _offsets[idx + 1]);
	}

	span<Vec> operator[](uint idx) {
		return span<Vec>(_vertices.data() + _offsets[idx], _vertices.data() + _offsets[idx + 1]);
	}

	uint size() const { return _offsets.size() - 1; }

	void reserve(uint rings, uint vertices) {
		_vertices.reserve(vertices);
		_offsets.reserve(rings + 1);
	}

	span<const Vec> vertices() const { return span<const Vec>(_vertices); }

	span<Vec> vertices() { return span<Vec>(_vertices); }

private:
	aligned_vector<Vec> _vertices;
	vector<uint> _offsets;
};

using xpolygon2 = xpolygon<double2>;
using xpolygon3 = xpolygon<double4>;

inline void format_e(string& s, string_view spec, const xpolygon2& p) {
	s += "(";
	for (size_t i = 0; i < p.size(); i++) {
		format_e(s, spec, p[i]);
		if (i != 0)
			s += " ";
	}
	s += ")";
}

inline void format_e(string& s, string_view spec, const xpolygon3& p) {
	s += "(";
	for (size_t i = 0; i < p.size(); i++) {
		format_e(s, spec, p[i]);
		if (i != 0)
			s += " ";
	}
	s += ")";
}

inline array<span<const double2>, 1> Rings(const polygon2& poly) {
	return { poly };
}

struct rings_iter {
	int ring;
	const xpolygon2& poly;

	optional<span<const double2>> next() {
		if (ring == poly.size())
			return nullopt;
		return poly[ring++];
	}
};

inline auto Rings(const xpolygon2& poly) {
	return iterable(rings_iter{ 0, poly });
}

template<typename Vec>
struct xpolygon_edge_iter {
	int ring = 0;
	int vertex = 0;
	const xpolygon<Vec>& poly;

	xpolygon_edge_iter(const xpolygon<Vec>& poly) : poly(poly) { }

	optional<segment<Vec>> next() {
		if (vertex >= poly[ring].size()) {
			vertex = 0;
			ring += 1;
		}
		if (ring >= poly.size())
			return {};
		span<const Vec> r = poly[ring];
		ON_SCOPE_EXIT(vertex += 1);
		return (vertex == 0) ? segment(r.back(), r[0]) : segment(r[vertex - 1], r[vertex]);
	}
};

template<typename T>
constexpr auto Edges(const xpolygon<T>& poly) {
	return iterable(xpolygon_edge_iter<T>(poly));
}

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

// Triangle meshes
using mesh2 = vector<triangle2>;
using mesh3 = aligned_vector<triangle3>;

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

inline double signed_double_area(const polygon2& poly) {
	double area = 0;
	if (poly.size() > 0)
		for (auto [a, b] : Edges(poly))
			area += signed_double_edge_area(a, b);
	return area;
}

inline double area(const polygon2& poly) {
	return abs(signed_double_area(poly));
}

template<typename Vec>
Vec AnyVertex(const vector<Vec>& p) { return p[0]; }
template<typename Vec>
Vec AnyVertex(const xpolygon<Vec>& p) { return p[0][0]; }
