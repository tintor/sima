#pragma once
#include <core/exception.h>
#include <geom/triangle.h>

// TODO move polygon and mesh stuff outside

using polygon2 = vector<double2>;
using polygon3 = vector<double3>;

// polygon2 with holes
template <typename Vec>
class xpolygon {
   public:
    xpolygon() { _offsets.push_back(0); }

    void add(Vec p) {
        assert(size() > 0);
        _vertices.push_back(p);
        _offsets.back() += 1;
    }

    void add(cspan<Vec> poly) {
        _vertices << poly;
        _offsets.push_back(_offsets.back() + poly.size());
    }

    cspan<Vec> operator[](uint idx) const {
        return {_vertices.data() + _offsets[idx], _offsets[idx + 1] - _offsets[idx]};
    }

    span<Vec> operator[](uint idx) { return {_vertices.data() + _offsets[idx], _offsets[idx + 1] - _offsets[idx]}; }

    uint size() const { return _offsets.size() - 1; }

    void reserve(uint rings, uint vertices) {
        _vertices.reserve(vertices);
        _offsets.reserve(rings + 1);
    }

    cspan<Vec> vertices() const { return cspan<Vec>(_vertices); }

    span<Vec> vertices() { return span<Vec>(_vertices); }

   private:
    vector<Vec> _vertices;
    vector<uint> _offsets;
};

using xpolygon2 = xpolygon<double2>;
using xpolygon3 = xpolygon<double3>;

inline void format_e(string& s, string_view spec, const xpolygon2& p) {
    s += "(";
    for (size_t i = 0; i < p.size(); i++) {
        format_e(s, spec, p[i]);
        if (i != 0) s += " ";
    }
    s += ")";
}

inline void format_e(string& s, string_view spec, const xpolygon3& p) {
    s += "(";
    for (size_t i = 0; i < p.size(); i++) {
        format_e(s, spec, p[i]);
        if (i != 0) s += " ";
    }
    s += ")";
}

inline array<cspan<double2>, 1> Rings(const polygon2& poly) { return {poly}; }

struct rings_iter {
    int ring;
    const xpolygon2& poly;

    optional<cspan<double2>> next() {
        if (ring == poly.size()) return nullopt;
        return poly[ring++];
    }
};

inline auto Rings(const xpolygon2& poly) { return iterable(rings_iter{0, poly}); }

template <typename Vec>
struct xpolygon_edge_iter {
    int ring = 0;
    int vertex = 0;
    const xpolygon<Vec>& poly;

    xpolygon_edge_iter(const xpolygon<Vec>& poly) : poly(poly) {}

    optional<segment<Vec>> next() {
        if (vertex >= poly[ring].size()) {
            vertex = 0;
            ring += 1;
        }
        if (ring >= poly.size()) return {};
        cspan<Vec> r = poly[ring];
        ON_SCOPE_EXIT(vertex += 1);
        return (vertex == 0) ? segment(r.back(), r[0]) : segment(r[vertex - 1], r[vertex]);
    }
};

template <typename T>
constexpr auto Edges(const xpolygon<T>& poly) {
    return iterable(xpolygon_edge_iter<T>(poly));
}

// TODO remove
template <typename T>
inline bool equal(const vector<T>& a, const vector<T>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (!equal(a[i], b[i])) return false;
    return true;
}

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

inline double signed_double_area(cspan<double2> poly) {
    double area = 0;
    if (poly.size() >= 3) {
        area += signed_double_edge_area(poly.back(), poly[0]);
        for (uint i = 1; i < poly.size(); i++) area += signed_double_edge_area(poly[i - 1], poly[i]);
    }
    return area;
}

inline double area(cspan<double2> poly) { return abs(signed_double_area(poly) / 2); }

// assumes that holes are oriented opposite from non-holes
inline double area(const xpolygon2& poly) {
    double sda = 0;
    for (auto ring : Rings(poly)) sda += signed_double_area(ring);
    return abs(sda / 2);
}

template <typename Vec>
Vec AnyVertex(const vector<Vec>& p) {
    return p[0];
}
template <typename Vec>
Vec AnyVertex(const xpolygon<Vec>& p) {
    return p[0][0];
}
