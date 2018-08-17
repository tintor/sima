#include "aabb.h"
#include "util.h"

/*pair<int4, int4> compute_aabb(span<const int4> vectors) {
	if (vectors.size() <= 1) {
		if (vectors.size() == 1)
			return {vectors[0], vectors[0]};
		int min = std::numeric_limits<int>::max();
		int max = std::numeric_limits<int>::min();
		return {int4{min, min, min, min}, int4{max, max, max, max}};
	}

	static_assert(sizeof(int4) * 2 == sizeof(int8));
	const int8* it;
	const int8* end;
	int8 min8, max8;

	assert(vectors.size() >= 2);
	if (aligned<32>(vectors.begin())) {
		it = reinterpret_cast<const int8*>(vectors.begin());
		end = it + vectors.size() / 2;
		min8 = *it++;
		max8 = min8;
	} else {
		assert(aligned<16>(vectors.begin()));
		it = reinterpret_cast<const int8*>(vectors.begin() + 1);
		end = it + (vectors.size() - 1) / 2;
		int8 v = vectors[0].xyzwxyzw;
		min8 = v;
		max8 = v;
	}
	assert(aligned<32>(it));
	assert(aligned<32>(end));

	while (it < end) {
		int8 v = *it;
		min8 = vmin(min8, v);
		max8 = vmax(max8, v);
		it += 1;
	}

	if (!aligned<32>(vectors.end())) {
		int8 v = vectors.back().xyzwxyzw;
		min8 = vmin(min8, v);
		max8 = vmax(max8, v);
	}

	int4 min4 = vmin(min8.xyzw, vshuffle(min8, min8, 4, 5, 6, 7));
	int4 max4 = vmax(max8.xyzw, vshuffle(max8, max8, 4, 5, 6, 7));
	return {min4, max4};
}*/
