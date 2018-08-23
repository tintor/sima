#pragma once
#include "span.h"
#include "vector.h"
#include "align_alloc.h"

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

using face = span<const point3>;

class xmesh3 {
public:
	xmesh3() { _offsets.push_back(0); }

	void add(face f) {
		_offsets.push_back(_offsets.back() + f.size());
		_vertices.resize(_vertices.size() + f.size());
		std::copy(f.begin(), f.end(), _vertices.end() - f.size());
	}

	uint size() const { return _offsets.size() - 1; }

	span<const point3> operator[](uint idx) const {
		return span<const point3>(_vertices.data() + _offsets[idx], _vertices.data() + _offsets[idx + 1]);
	}

	// might not be needed for convex hull after all!
	void erase(uint idx) {
		uint s = _offsets[idx + 1] - _offsets[idx];
		for (uint i = _offsets[idx]; i < _vertices.size() - s; i++)
			_vertices[i] = _vertices[i + s];
		_vertices.resize(_vertices.size() - s);
		for (uint i = idx + 1; i < _offsets.size() - 1; i++)
			_offsets[i] = _offsets[i + 1] - s;
		_offsets.resize(_offsets.size() - 1);
	}

	class iterator {
		uint _pos;
		const xmesh3& _mesh;
	public:
		iterator(uint pos, const xmesh3& mesh) : _pos(pos), _mesh(mesh) { }
		face operator*() const { return _mesh[_pos]; }
		iterator& operator++() { _pos++; return *this; }
		bool operator!=(iterator other) const { return _pos != other._pos; }
	};

	iterator begin() const { return iterator(0, *this); }
	iterator end() const { return iterator(size(), *this); }

	span<const point3> vertices() const { return _vertices; }
	span<const uint> offsets() const { return _offsets; }

private:
	aligned_vector<point3> _vertices;
	vector<uint> _offsets;
};
