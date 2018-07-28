#pragma once
#include "range.h"
#include "triangle.h"
#include "array_ptr.h"

template<typename Vec>
struct aabb {
	Vec min, max;

	aabb() {
		reset();
	}

	auto range_dim() const { return range(sizeof(min) / sizeof(min.x)); }

	explicit aabb(Vec p) {
		min = max = p;
	}

	explicit aabb(Vec a, Vec b) {
		min = max = a;
		add(b);
	}

	explicit aabb(Vec a, Vec b, Vec c) {
		min = max = a;
		add(b);
		add(c);
	}

	explicit aabb(const segment<Vec>& p) {
		min = max = p.a;
		add(p.b);
	}

	explicit aabb(const triangle<Vec>& p) {
		min = max = p.a;
		add(p.b);
		add(p.c);
	}

	explicit aabb(array_cptr<Vec> p) {
		reset();
		for (auto v : p)
			add(v);
	}

	explicit aabb(array_cptr<triangle<Vec>> ap) {
		reset();
		for (auto p : ap) {
			add(p.a);
			add(p.b);
			add(p.c);
		}
	}

	bool operator==(const aabb<Vec>& v) const {
		return min == v.min && max == v.max;
	}

	bool operator!=(const aabb<Vec>& v) const {
		return !operator==(v);
	}

	void reset() {
		for (auto i : range_dim()) {
			min[i] = std::numeric_limits<decltype(min.x)>::max();
			max[i] = std::numeric_limits<decltype(min.x)>::min();
		}
	}

	void add(Vec v) {
		for (auto i : range_dim()) {
			if (v[i] < min[i])
				min[i] = v[i];
			if (v[i] > max[i])
				max[i] = v[i];
		}
	}

	bool valid() const {
		for (auto i : range_dim())
			if (min[i] > max[i])
				return false;
		return true;
	}

	Vec size() const {
		assert(valid());
		return max - min;
	}

	Vec center() const {
		assert(valid());
		return (min + max) / 2;
	}

	// strictly inside!
	bool inside(ivec3 e) const {
		for (auto i : range_dim())
			if (e[i] <= min[i] || e[i] >= max[i])
				return false;
		return true;
	}

	bool intersects(ivec3 e) const {
		for (auto i : range_dim())
			if (e[i] < min[i] || e[i] > max[i])
				return false;
		return true;
	}

	template<typename T>
	static bool overlaps(T amin, T amax, T bmin, T bmax) {
		return bmin < amax && amin < bmax;
	}

	template<typename T>
	static bool intersects(T amin, T amax, T bmin, T bmax) {
		return bmin <= amax && amin <= bmax;
	}

	bool overlaps(const aabb& e) const {
		for (auto i : range_dim())
			if (!overlaps(min[i], max[i], e.min[i], e.max[i]))
				return false;
		return true;
	}

	bool intersects(const aabb& e) const {
		for (auto i : range_dim())
			if (!intersects(min[i], max[i], e.min[i], e.max[i]))
				return false;
		return true;
	}
};

template<typename T>
void format_e(std::string& s, std::string_view spec, const aabb<T>& box) {
	s += "aabb(";
	format_e(s, "", box.min);
	s += ", ";
	format_e(s, "", box.max);
	s += ')';
}
