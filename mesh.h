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

template<typename T>
class array_iterator {
	uint _pos;
	const T& _array;
public:
	array_iterator(uint pos, const T& a) : _pos(pos), _array(a) { }
	auto operator*() const { return _array[_pos]; }
	array_iterator& operator++() { _pos++; return *this; }
	bool operator!=(array_iterator other) const { return _pos != other._pos; }
};

// face can have arbitrary number of vertices and holes
// first ring is exterior, other rings are interior and oriented opposite
struct face {
	face(const point3* vertices, span<const uint> offsets) : _vertices(vertices), _offsets(offsets) { }
	span<const point3> operator[](uint idx) const {
		return span<const point3>(_vertices + _offsets[idx], _vertices + _offsets[idx + 1]);
	}
	uint size() const { return _offsets.size() - 1; }
	uint vertex_size() const { return _offsets.back() - _offsets[0]; }

	span<const uint> offsets() const { return _offsets; }
	span<const point3> vertices() const {
		return span<const point3>(_vertices + _offsets[0], _vertices + _offsets.back());
	}

	auto begin() const { return array_iterator<face>(0, *this); }
	auto end() const { return array_iterator<face>(size(), *this); }
private:
	const point3* _vertices;
	span<const uint> _offsets;
};

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
		for (uint offset : f.offsets().pop_front())
			_offsets.push_back(offset + a);

		_ring_offsets.push_back(_ring_offsets.back() + f.size());
	}

	uint size() const { return _ring_offsets.size() - 1; }

	face operator[](uint idx) const {
		span<const uint> offsets(_offsets.data() + _ring_offsets[idx], _offsets.data() + _ring_offsets[idx + 1]);
		return face(_vertices.data(), offsets);
	}

	auto begin() const { return array_iterator<xmesh3>(0, *this); }
	auto end() const { return array_iterator<xmesh3>(size(), *this); }

	span<const point3> vertices() const { return _vertices; }

private:
	aligned_vector<point3> _vertices;
	vector<uint> _offsets;
	vector<uint> _ring_offsets;
};

// +1 - disjoint
//  0 - touching / contact
// -1 - overlaping / penetrating
int Classify(const xmesh3& ma, const xmesh3& mb);
