#include <geom/broadphase.h>

void HashBroadphase::getCandidates(sphere s, vector<int>& candidates) const {
	if (_grid.empty())
		return;
	double4 r = broad4(s.radius() + _max_radius);
	double4 a = floor((s.center() - r) * _inv_cell);
	double4 b = ceil((s.center() + r) * _inv_cell);
	int ax = a.x, ay = a.y, az = a.z;
	int bx = b.x, by = b.y, bz = b.z;

	for (int x = ax; x <= bx; x++)
		for (int y = ay; y <= by; y++)
			for (int z = az; z <= bz; z++) {
				ulong cell = encode(x, y, z);
				auto b = _grid.bucket(cell);
				auto end = _grid.end(b);
				for (auto it = _grid.begin(b); it != end; it++)
					candidates.push_back(it->second);
			}
}
