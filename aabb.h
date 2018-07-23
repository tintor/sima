#pragma once
#include "range.h"
#include "triangle.h"
#include "array_ptr.h"

template<typename T>
struct interval {
	T min = std::numeric_limits<T>::max();
	T max = std::numeric_limits<T>::min();

	void include(T a) {
		if (a < min)
			min = a;
		if (a > max)
			max = a;
	}

	bool overlaps(interval e) const {
		return e.min < max && min < e.max;
	}
	
	bool intersects(interval e) const {
		return e.min <= max && min <= e.max;
	}
};

template<typename Vec>
struct aabb {
	static constexpr int dim = sizeof(Vec) / sizeof(Vec::x);
	std::array<interval<decltype(Vec::x)>, dim> mm;
	
	explicit aabb(Vec p) {
		for (auto i : range(dim))
			mm[i].include(p[i]);
	}

	explicit aabb(Vec a, Vec b) {
		for (auto i : range(dim)) {
			mm[i].include(a[i]);
			mm[i].include(b[i]);
		}
	}

	explicit aabb(Vec a, Vec b, Vec c) {
		for (auto i : range(dim)) {
			mm[i].include(a[i]);
			mm[i].include(b[i]);
			mm[i].include(c[i]);
		}
	}

	explicit aabb(const segment<Vec>& p) {
		for (auto i : range(dim)) {
			mm[i].include(p.a[i]);
			mm[i].include(p.b[i]);
		}
	}

	explicit aabb(const triangle<Vec>& p) {
		for (auto i : range(dim)) {
			mm[i].include(p.a[i]);
			mm[i].include(p.b[i]);
			mm[i].include(p.c[i]);
		}
	}

	explicit aabb(array_cptr<Vec> p) {
		for (auto v : p)
			for (auto i : range(dim))
				mm[i].include(v[i]);
	}

	explicit aabb(array_cptr<triangle<Vec>> ap) {
		for (auto i : range(dim))
			for (auto p : ap) {
				mm[i].include(p.a[i]);
				mm[i].include(p.b[i]);
				mm[i].include(p.c[i]);
			}
	}

	// strictly inside!
	bool inside(ivec3 e) const {
		for (auto i : range(dim))
			if (e[i] <= mm[i].min || e[i] >= mm[i].max)
				return false;
		return true;
	}

	bool intersects(ivec3 e) const {
		for (auto i : range(dim))
			if (e[i] < mm[i].min || e[i] > mm[i].max)
				return false;
		return true;
	}

	bool overlaps(const aabb& e) const {
		for (auto i : range(dim))
			if (!mm[i].overlaps(e.mm[i]))
				return false;
		return true;
	}

	bool intersects(const aabb& e) const {
		for (auto i : range(dim))
			if (!mm[i].intersects(e.mm[i]))
				return false;
		return true;
	}
};
