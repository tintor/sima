#pragma once
#include "mesh.h"
#include "segment.h"
#include "plane.h"
#include "aabb.h"

inline aabb2 Box(const polygon2& p) { return aabb2(p).buffer(Tolerance); }
inline aabb2 Box(const xpolygon2& p) { return aabb2(p.vertices()).buffer(Tolerance); }

// All Sign() functions return:
// +1 above / positive side
//  0 contact / between
// -1 below / negative side

inline int Sign(double d) {
	if (d > Tolerance) return 1;
	if (d < -Tolerance) return -1;
	return 0;
}

// TESTED
inline int Sign(segment2 s, double2 v) {
	return Sign(signed_double_area(s.a, s.b, v) / length(s.a - s.b));
}

inline int Sign(ray2 s, double2 v) {
	return Sign(signed_double_area(s.origin, s.origin + s.unit_dir, v));
}

inline int Sign(plane p, double4 v) { return Sign(p.distance(v)); }

// All Classify() functions return:
// +1 disjoint / separate
//  0 touching
// -1 overlap / penetrating

// TODO move templates inside classify.cc for faster compilation
template<typename Polygon2>
int Classify(const Polygon2& a, double2 p) {
	return Classify(a, p, Box(a));
}

template<typename Polygon2>
int Classify(const Polygon2& a, double2 p, aabb2 box) {
	if (!box.intersects(p))
		return +1;

	int crossings = 0; // 1 - half a crossing, 2 - full crossing
	ray2 ray(p, p + double2{1, 0});
	for (auto e : Edges(a)) {
		int say = Sign(e.a.y - p.y);
		int sby = Sign(e.b.y - p.y);
		if (say * sby > 0) // edge is above or below ray
			continue;
		if (say == 0 && sby == 0) { // colinear
			int sax = Sign(e.a.x - p.x);
			int sbx = Sign(e.b.x - p.x);
			if (sax * sbx <= 0) // TODO this will make vertices square instead of round
				return 0;
		}
		if (sby == 0) { // in edge
			crossings += say;
			continue;
		}
		if (say == 0) { // out edge
			crossings -= sby;
			continue;
		}
		// E is not horizontal
		int sp = Sign(e, p); // TODO check proper sign!
		if (sp == 0)
			return 0;
		if (sp < 0) // P is on the left side of e
			crossings += 2;
	}
	return ((crossings / 2) & 1) ? -1 : 1;
}

class Intervals {
public:
	void add(double begin, double end) {
		_points.emplace_back(true, begin);
		_points.emplace_back(false, end);
	}
	void unionAll();
	size_t size() const { return _points.size() / 2; }
	pair<double, double> operator[](size_t idx) {
		return { _points[idx * 2].second, _points[idx * 2 + 1].second };
	}
private:
	vector<pair<bool, double>> _points;
};

template<typename Polygon2>
int Classify(const Polygon2& f, segment2 s, vector<pair<double, double>>* intersections = nullptr) {
	return Classify(f, s, Box(f), intersections);
}

template<typename Polygon2>
int Classify(const Polygon2& f, segment2 s, aabb2 box, vector<pair<double, double>>* intersections = nullptr) {
	if (!Intersects(box, aabb2(s)))
		return +1;

	if (!intersections)
		if (Classify(f, s.a, box) == -1 || Classify(f, s.b, box) == -1)
			return -1;

	Intervals intervals;
	for (segment2 e : Edges(f)) {
		double2 t;
		if (relate(e, s, nullptr, &t) != 'D') {
			double x = t.x, y = t.y;
			if (intersections)
				intersections->emplace_back(x, y);
			intervals.add(x, y);
		}
	}

	if (Classify(f, s.a, box) == -1 || Classify(f, s.b, box) == -1)
		return -1;
	if (intervals.size() == 0)
		return +1;

	// check all mid points between travels if they are inside polygon
	intervals.unionAll();
	for (size_t i = 1; i < intervals.size(); i++) {
		double t = (intervals[i - 1].second + intervals[i].first) / 2;
		int c = Classify(f, s.linear(t), box);
		assert(c != 0);
		if (c == -1)
			return -1;
	}
	return 0;
}

template<typename Polygon2>
int Classify(const Polygon2& a, const Polygon2& b) {
	aabb2 va = Box(a), vb = Box(b);
	if (!Intersects(va, vb))
		return +1;

	if (Classify(a, AnyVertex(b), va) < 0 || Classify(b, AnyVertex(a), vb) < 0)
		return -1;

	int result = 1;
	for (auto ea : Edges(a)) {
		int c = Classify(b, ea, vb);
		if (c == -1)
			return -1;
		if (c == 0)
			result = 0;
	}
	for (auto eb : Edges(b)) {
		int c = Classify(a, eb, va);
		if (c == -1)
			return -1;
		if (c == 0)
			result = 0;
	}
	return result;
}

int Classify(const xpolygon2& f, ray2 s);

int Classify(const face& f, double4 v);
int Classify(const face& f, const ray3& s, double* travel = nullptr);
pair<int, int> ClassifyDoubleSided(const face& f, const ray3& s);
int Classify(plane p, const segment3& s);
int Classify(const face& f, const segment3& s, vector<pair<double, double>>* intersections = nullptr);
int Classify(const xmesh3& m, double4 p, const aabb4& box);
int Classify(const xmesh3& m, const segment3& s, const aabb4& box);
int Classify(const xmesh3& ma, const xmesh3& mb);
