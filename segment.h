#pragma once
#include "vector.h"
#include "hash.h"
#include <cassert>

// finite part of line between two points
template<typename Vec>
struct segment {
	Vec a, b;

	constexpr segment(Vec a, Vec b) : a(a), b(b) { }
	constexpr segment(const segment& s) : a(s.a), b(s.b) { }

	constexpr bool operator==(segment v) const { return equal(a, v.a) && equal(b, v.b); }
	constexpr bool operator!=(segment v) const { return !(operator==(v)); }

	segment reversed() const { return {b, a}; }

	Vec linear(double t) const { return a + (b - a) * t; }

	// Inverse of linear
	double param(Vec p) const {
		Vec d = b - a;
		return dot(p - a, d) / dot(d, d);
	}

	// Return point on segment that is nearest to P
	Vec nearest(Vec p) const {
		Vec d = b - a;
		double t = dot(p - a, d);
		if (t <= 0)
			return a;
		if (t >= dot(d, d))
			return b;
		return a + d * (t / dot(d, d));
	}
};

using segment2 = segment<double2>;
using segment3 = segment<point3>;

// starts from a, and goes to infinity
struct ray3 {
	const double4 unit_dir;
   	const point3 origin;

	ray3(const ray3& s) : unit_dir(s.unit_dir), origin(s.origin) { }
	ray3(segment3 s) : unit_dir(normalize(s.b - s.a)), origin(s.a) { }
	ray3(double4 a, double4 b) : unit_dir(normalize(b - a)), origin(a) { }

	bool operator==(ray3 v) const { return equal(origin, v.origin) && equal(unit_dir, v.unit_dir); }
	bool operator!=(ray3 v) const { return !(operator==(v)); }

	double4 linear(double t) const { return origin + unit_dir * t; }

	// Inverse of linear
	double param(double4 p) const {
		return dot(p - origin, unit_dir);
	}
};

// TODO raw_line3: which just stores points A and B (no normalization)

// infinite on both sides unlike segment
// Note: dir may be flipped after construction!
struct line3 {
	const double4 unit_dir, origin;

	static double4 normalize_dir(double4 d) {
		if (d.x < 0)
			return -d;
		if (d.x == 0 && d.y < 0)
			return -d;
		if (d.x == 0 && d.y == 0 && d.z < 0)
			return -d;
		return d;
	}

	line3(const line3& s) : unit_dir(s.unit_dir), origin(s.origin) { }
	line3(ray3 s)
		: unit_dir(normalize_dir(s.unit_dir))
		, origin(s.origin - unit_dir * dot(s.origin, unit_dir)) {
	}
	line3(segment3 s) : line3(ray3(s)) { }
	line3(double4 a, double4 b) : line3(segment3{a, b}) { }

	// TODO(marko) equality condition should check if lines are overlapping
	bool operator==(line3 v) const { return equal(origin, v.origin) && equal(unit_dir, v.unit_dir); }
	bool operator!=(line3 v) const { return !(operator==(v)); }

	double4 linear(double t) const { return origin + unit_dir * t; }

	double4 nearest(double4 p) const {
		return origin + unit_dir * dot(p - origin, unit_dir);
	}
};

// Tested and stable!
enum class NearestCase {
	RayRay = 0,
	RayPoint = 1,
	PointPoint = 2,
};

pair<segment3, NearestCase> nearest(segment3 p, segment3 q);

inline void format_e(string& s, string_view spec, line3 v) {
	format_s(s, "line(%s, %s)", v.origin, v.unit_dir);
}

inline void format_e(string& s, string_view spec, ray3 v) {
	format_s(s, "ray(%s, %s)", v.origin, v.unit_dir);
}

inline void format_e(string& s, string_view spec, segment3 v) {
	format_s(s, "segment(%s, %s)", v.a, v.b);
}

namespace std {

template<> struct hash<line3> {
	size_t operator()(line3 s) const { return ::hash(s.origin, s.unit_dir); }
};

template<> struct hash<ray3> {
	size_t operator()(ray3 s) const { return ::hash(s.origin, s.unit_dir); }
};

template<typename Vec> struct hash<segment<Vec>> {
	size_t operator()(segment<Vec> s) const { return ::hash(s.a, s.b); }
};

}

inline auto angle(segment3 p, segment3 q) {
	return angle(p.b - p.a, q.b - q.a);
}
