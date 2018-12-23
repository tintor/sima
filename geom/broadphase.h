#include <core/std.h>
#include <geom/sphere.h>

class HashBroadphase {
public:
	void reset(double cell_size, double max_radius) {
		_inv_cell = 1 / cell_size;
		_max_radius = max_radius;
		clear();
	}

	void getCandidates(sphere s, vector<int>& candidates) const;

	void add(sphere s, int id) {
		double4 c = floor(s.center() * _inv_cell);
		int x = c.x;
		int y = c.y;
		int z = c.z;
		ulong cell = encode(x, y, z);
		_grid.insert({cell, id});
	}

	void clear() { _grid.clear(); }
private:
	static ulong encode(ulong x, ulong y, ulong z) {
		return ((x & 0x1FFFFF) << 42) | ((y & 0x1FFFFF) << 21) | (z & 0x1FFFFF);
	}

	double _inv_cell;
	double _max_radius;
	unordered_multimap<ulong, int> _grid;
};
