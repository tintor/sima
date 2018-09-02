#pragma once
#include "vector.h"
#include "hash.h"
#include <cassert>

constexpr double Tolerance = 0.5e-6;

template<typename Vec>
bool Equals(Vec a, Vec b) { return squared(a - b) <= squared(Tolerance); }

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
using segment3 = segment<double4>;

// starts from a, and goes to infinity
template<typename Vec>
struct ray {
	const Vec unit_dir, origin;

	constexpr ray(const ray& s) : unit_dir(s.unit_dir), origin(s.origin) { }
	ray(segment<Vec> s) : unit_dir(normalize(s.b - s.a)), origin(s.a) { }
	ray(Vec a, Vec b) : unit_dir(normalize(b - a)), origin(a) { }

	// using max() instead of infinity() in order for mult with zero to stay zero
	Vec infinity() { return linear(std::numeric_limits<double>::max()); }

	bool operator==(ray v) const { return equal(origin, v.origin) && equal(unit_dir, v.unit_dir); }
	bool operator!=(ray v) const { return !(operator==(v)); }

	Vec linear(double t) const { return origin + unit_dir * t; }

	// Inverse of linear
	double param(Vec p) const {
		return dot(p - origin, unit_dir);
	}
};

using ray2 = ray<double2>;
using ray3 = ray<double4>;

// TODO raw_line3: which just stores points A and B (no normalization)

// infinite on both sides unlike segment
// Note: dir may be flipped after construction!
template<typename Vec>
struct line {
	const Vec unit_dir, origin;

	static Vec normalize_dir(Vec d) {
		if (d.x < 0)
			return -d;
		if (d.x == 0 && d.y < 0)
			return -d;
		if (d.x == 0 && d.y == 0 && d.z < 0)
			return -d;
		return d;
	}

	line(const line& s) : unit_dir(s.unit_dir), origin(s.origin) { }
	line(ray<Vec> s)
		: unit_dir(normalize_dir(s.unit_dir))
		, origin(s.origin - unit_dir * dot(s.origin, unit_dir)) {
	}
	line(segment<Vec> s) : line(ray(s)) { }
	line(Vec a, Vec b) : line(segment<Vec>{a, b}) { }

	// TODO(marko) equality condition should check if lines are overlapping
	bool operator==(line<Vec> v) const { return equal(origin, v.origin) && equal(unit_dir, v.unit_dir); }
	bool operator!=(line<Vec> v) const { return !(operator==(v)); }

	Vec linear(double t) const { return origin + unit_dir * t; }

	Vec nearest(Vec p) const {
		return origin + unit_dir * dot(p - origin, unit_dir);
	}
};

using line2 = line<double2>;
using line3 = line<double4>;

char relate(segment<double2> p, segment<double2> q);
bool relate(segment<double2> p, ray<double2> q);

// Tested and stable!
enum class NearestCase {
	RayRay = 0,
	RayPoint = 1,
	PointPoint = 2,
};

pair<segment3, NearestCase> nearest(segment3 p, segment3 q);

template<typename Vec>
void format_e(string& s, string_view spec, line<Vec> v) {
	format_s(s, "line(%s, %s)", v.origin, v.unit_dir);
}

template<typename Vec>
void format_e(string& s, string_view spec, ray<Vec> v) {
	format_s(s, "ray(%s, %s)", v.origin, v.unit_dir);
}

template<typename Vec>
void format_e(string& s, string_view spec, segment<Vec> v) {
	format_s(s, "segment(%s, %s)", v.a, v.b);
}

template<typename Vec>
hash operator<<(hash h, line<Vec> s) { return h << s.origin << s.unit_dir; }
template<typename Vec>
hash operator<<(hash h, ray<Vec> s) { return h << s.origin << s.unit_dir; }
template<typename Vec>
hash operator<<(hash h, segment<Vec> s) { return h << s.a << s.b; }

inline auto angle(segment3 p, segment3 q) {
	return angle(p.b - p.a, q.b - q.a);
}

inline double signed_double_edge_area(double2 a, double2 b) {
	return (a.x + b.x) * (a.y - b.y);
}

inline double signed_double_area(double2 a, double2 b, double2 c) {
	return signed_double_edge_area(a, b) + signed_double_edge_area(b, c) + signed_double_edge_area(c, a);
}

inline int Classify(double d) {
	if (d > Tolerance) return 1;
	if (d < -Tolerance) return -1;
	return 0;
}

inline int Classify(segment2 s, double2 v) {
	return Classify(signed_double_area(s.a, s.b, v) / length(s.a - s.b));
}

inline int Classify(ray2 s, double2 v) {
	return Classify(signed_double_area(s.origin, s.origin + s.unit_dir, v));
}

inline bool Colinear(line3 s, double4 v) {
	auto va = v - s.origin;
    return squared(va - s.unit_dir * dot(va, s.unit_dir)) <= squared(Tolerance);
}

char relate(segment2 p, segment2 q);

inline bool relate_abxo(segment2 p, segment2 q) {
	char c = relate(p, q);
	return c == 'A' || c == 'B' || c == 'X' || c == 'O';
}
