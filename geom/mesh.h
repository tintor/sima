#pragma once
#include <core/auto.h>
#include <core/each.h>
#include <core/exception.h>
#include <core/std.h>
#include <geom/plane.h>
#include <geom/segment.h>
#include <geom/vector.h>

using mesh2 = vector<triangle2>;
using mesh3 = vector<triangle3>;

inline string wkt(const mesh2& mesh) {
    string s;
    s += "MULTIPOLYGON ((";
    for (auto i : range(mesh.size())) {
        if (i != 0) s += ", ";
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

// TODO should we store plane for every face?
// store normalized plane or unnormalized plane?
// pros:
// - faster point_in_convex
// - faster distance_from_convex
// cons:
// - slower transform

// used when loading meshes from file
//
// used when constructing meshes from convex hull
//
// used when constructing meshes directly

template <typename T>
class array_iterator {
    uint _pos;
    const T& _array;

   public:
    array_iterator(uint pos, const T& a) : _pos(pos), _array(a) {}
    auto operator*() const { return _array[_pos]; }
    array_iterator& operator++() {
        _pos++;
        return *this;
    }
    bool operator!=(array_iterator other) const { return _pos != other._pos; }
};

// face can have arbitrary number of vertices and holes
// first ring is exterior, other rings are interior and oriented opposite
struct face {
    face(const double3* vertices, cspan<uint> offsets, plane p) : _vertices(vertices), _offsets(offsets), _plane(p) {}
    cspan<double3> operator[](uint idx) const {
        return {_vertices + _offsets[idx], _offsets[idx + 1] - _offsets[idx]};
    }
    plane plane() const { return _plane; }
    uint size() const { return _offsets.size() - 1; }
    uint vertex_size() const { return _offsets.back() - _offsets[0]; }

    cspan<uint> offsets() const { return _offsets; }
    cspan<double3> vertices() const { return {_vertices + _offsets[0], _offsets.back() - _offsets[0]}; }

    auto begin() const { return array_iterator<face>(0, *this); }
    auto end() const { return array_iterator<face>(size(), *this); }

   private:
    const double3* _vertices;
    cspan<uint> _offsets;
    class plane _plane;
};

struct face_edge_iter {
    int ring = 0;
    int vertex = 0;
    const face& f;

    face_edge_iter(const face& f) : f(f) {}

    optional<segment3> next() {
        if (vertex >= f[ring].size()) {
            vertex = 0;
            ring += 1;
        }
        if (ring >= f.size()) return {};
        cspan<double3> r = f[ring];
        ON_SCOPE_EXIT(vertex += 1);
        return (vertex == 0) ? segment(r.back(), r[0]) : segment(r[vertex - 1], r[vertex]);
    }
};

inline auto Edges(const face& f) { return iterable(face_edge_iter(f)); }

class xmesh3 {
   public:
    xmesh3() {
        _offsets.push_back(0);
        _ring_offsets.push_back(0);
    }

    void add(face f) {
        _vertices << f.vertices();

        _offsets.reserve(_offsets.size() + f.offsets().size() - 1);
        uint a = _offsets.back();
        for (uint offset : f.offsets().subspan(1)) _offsets.push_back(offset + a);

        _ring_offsets.push_back(_ring_offsets.back() + f.size());
    }

    uint size() const { return _ring_offsets.size() - 1; }

    face operator[](uint idx) const {
        cspan<uint> offsets(_offsets.data() + _ring_offsets[idx], _ring_offsets[idx + 1] - _ring_offsets[idx]);
        return face(_vertices.data(), offsets, _facePlanes[idx]);
    }

    auto begin() const { return array_iterator<xmesh3>(0, *this); }
    auto end() const { return array_iterator<xmesh3>(size(), *this); }

    // contains duplicates!
    cspan<double3> vertices() const { return _vertices; }

    cspan<double3> uniqueVertices() const { THROW(not_implemented); }
    vector<pair<segment3, double3>> uniqueEdges() const { THROW(not_implemented); }

   private:
    vector<double3> _vertices;
    vector<uint> _offsets;
    vector<uint> _ring_offsets;

    vector<plane> _facePlanes;
    vector<segment3> _uniqueEdges;
    vector<double3> _edgeNormals;
    vector<double3> _uniqueVertices;
};
