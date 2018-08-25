#pragma once
#include "range.h"
#include "triangle.h"
#include "span.h"

pair<int4, int4> compute_aabb(span<const int4> vectors);

template<typename Vec>
struct aabb {
	Vec min, max;

	aabb() {
		reset();
	}

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

	explicit aabb(segment<Vec> p) {
		min = max = p.a;
		add(p.b);
	}

	explicit aabb(triangle<Vec> p) {
		min = max = p.a;
		add(p.b);
		add(p.c);
	}

	explicit aabb(span<const Vec> p) {
		reset();
		for (auto v : p)
			add(v);
	}

	explicit aabb(span<const triangle<Vec>> ap) {
		reset();
		for (auto p : ap) {
			add(p.a);
			add(p.b);
			add(p.c);
		}
	}

	bool operator==(aabb v) const { return equal(min, v.min) && equal(max, v.max); }
	bool operator!=(aabb v) const { return !operator==(v); }

	void reset() {
		broadcast(min, std::numeric_limits<double>::max());
		broadcast(max, std::numeric_limits<double>::min());
	}

	void add(Vec v) {
		min = vmin(v, min);
		max = vmax(v, max);
	}

	bool valid() const {
		return all(min <= max);
	}

	auto size() const {
		assert(valid());
		return max - min;
	}

	Vec center() const {
		assert(valid());
		return (min + max) / 2;
	}

	bool intersects(Vec e) const {
		// TODO vector intrinsics
		for (auto i : range(3))
			if (e[i] < min[i] || e[i] > max[i])
				return false;
		return true;
	}

	template<typename T>
	static bool overlaps(T amin, T amax, T bmin, T bmax) {
		return bmin < amax && amin < bmax;
	}

	bool overlaps(aabb e) const {
		// TODO vector intrinsics
		for (auto i : range(3))
			if (!overlaps(min[i], max[i], e.min[i], e.max[i]))
				return false;
		return true;
	}

	template<typename T>
	static bool intersects(T amin, T amax, T bmin, T bmax) {
		return bmin <= amax && amin <= bmax;
	}

	bool intersects(aabb e) const {
		// TODO vector intrinsics
		for (auto i : range(3))
			if (!intersects(min[i], max[i], e.min[i], e.max[i]))
				return false;
		return true;
	}
};

// double only!
using aabb2 = aabb<double2>;
using aabb4 = aabb<double4>;

template<typename Vec>
void format_e(string& s, string_view spec, aabb<Vec> box) {
	s += "aabb(";
	format_e(s, spec, box.min);
	s += ", ";
	format_e(s, spec, box.max);
	s += ')';
}

// uniform inside a box (not very efficient)
template<typename RND>
double2 uniform2(RND& rnd, aabb2 box) {
	auto s = box.size();
	double m = max(s.x, s.y);
	while (true) {
		double2 v = uniform2(rnd, 0, m) + box.min;
		if (box.intersects(v))
			return v;
	}
}

// uniform inside a box (not very efficient)
template<typename RND>
double4 uniform3(RND& rnd, aabb4 box) {
	auto s = box.size();
	double m = max(s.x, s.y, s.z);
	while (true) {
		double4 v = uniform3v(rnd, 0, m) + box.min;
		if (box.intersects(v))
			return v;
	}
}
