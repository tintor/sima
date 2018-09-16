#include "classify.h"
#include "exception.h"
#include "edges.h"
#include "project.h"
#include "properties.h"

int Classify(const xpolygon2& f, ray2 s) {
	THROW(not_implemented);
}

pair<int, int> ClassifyDoubleSided(const xpolygon2& f, ray2 s) {
	THROW(not_implemented);
}

int Classify(const face& f, double4 v) {
	if (!Intersects(aabb4(f.vertices()), aabb4(v)))
		return +1;

	plane p = best_fit_plane(f.vertices());
	double d = p.distance(v);
	if (abs(d) > Tolerance)
		return +1;

	int axis = ProjectionAxis(f);
	int c = Classify(Project(f, axis), Project(v - p.normal() * d, axis));
	return (c == 1) ? 1 : 0;
}

int Classify(const face& f, const ray3& s, double* travel) {
	if (!Intersects(aabb4(f.vertices()), aabb4(s.origin, s.infinity())))
		return +1;

	plane p = best_fit_plane(f.vertices());

	double d = p.distance(s.origin);
	double t = -d / p.distance(s.unit_dir);
	if (travel)
		*travel = t;
	if (std::isfinite(t)) {
		if (t < -Tolerance)
			return +1;
		return Classify(f, s.linear(t));
	}

	if (abs(d) > Tolerance)
		return +1;

	int axis = ProjectionAxis(f);
	int c = Classify(Project(f, axis), Project(s, axis));
	return (c == 1) ? 1 : 0;
}

// Computes Classify for both s and -s at once.
pair<int, int> ClassifyDoubleSided(const face& f, const ray3& s) {
	// TODO bounding box
	//if (!Intersects(aabb4(f.vertices()), aabb4(s))
	//	return +1;

	plane p = best_fit_plane(f.vertices());

	double d = p.distance(s.origin);
	double t = -d / p.distance(s.unit_dir);
	if (std::isfinite(t)) {
		int c = Classify(f, s.linear(t));
		if (c == 1)
			return {1, 1};
		if (t > Tolerance)
			return {c, 1};
		if (t < -Tolerance)
			return {1, c};
		return {c, c};
	}

	if (d > Tolerance || d < -Tolerance)
		return {1, 1};

	int axis = ProjectionAxis(f);
	auto c = ClassifyDoubleSided(Project(f, axis), Project(s, axis));
	if (c.first == -1)
		c.first = 0;
	if (c.second == -1)
		c.second = 0;
	return c;
}

int Classify(plane p, const segment3& s) {
	int ca = Sign(p, s.a);
	int cb = Sign(p, s.b);
	if (ca == +1 && cb == +1)
		return +1;
	if (ca == -1 && cb == -1)
		return +1;
	if (ca == 0 || cb == 0)
		return 0;
	return -1;
}

// [intersections] represent relative travels along segment S
int Classify(const face& f, const segment3& s, vector<pair<double, double>>* intersections) {
	if (!Intersects(aabb4(f.vertices()), aabb4(s)))
		return +1;

	plane p = best_fit_plane(f.vertices());
	double da = p.distance(s.a);
	double db = p.distance(s.b);
	int ca = Sign(da);
	int cb = Sign(db);

	if (ca == +1 && cb == +1)
		return +1;
	if (ca == -1 && cb == -1)
		return +1;

	int axis = ProjectionAxis(f);
	auto f2 = Project(f, axis);

	if (ca == 0 && cb == 0) {
		// coplanar case
		auto s2 = Project(s, axis);
		return Classify(f2, s2, intersections);
	}

	if (ca == 0 || cb == 0) {
		// segment touches the plane
		auto v2 = Project((ca == 0) ? s.a : s.b, axis);
		if (intersections) {
			double t = (ca == 0) ? 0 : 1;
			intersections->emplace_back(t, t);
		}
		return std::max(0, Classify(f2, v2));
	}

	// segment intersects the plane
	double t = da / (da - db);
	auto v = s.linear(t);
	if (intersections) {
		intersections->emplace_back(t, t);
	}
	return Classify(f2, Project(v, axis));
}

int Classify(const xmesh3& m, double4 p, const aabb4& box) {
	if (!box.intersects(p))
		return +1;

	// Check if the point is on the boundary
	for (const face& f : m)
		if (Classify(f, p) == 0)
			return 0;

	// Shoot a ray (in both directions at once) and count crossings
	// (if both ray hits any edge or vertex, shoot a new ray and repeat)
	std::default_random_engine rnd;
	while (true) {
		ray3 ray(p, p + uniform_dir3(rnd));
		pair<int, int> result = {1, 1};
		for (const face& f : m) {
			auto c = ClassifyDoubleSided(f, ray);
			result.first *= c.first;
			result.second *= c.second;
			if (result.first == 0 && result.second == 0)
				break;
		}
		if (result.first == -1 || result.second == -1)
			return -1;
		if (result.first == 1 || result.second == 1)
			return 1;
	}
}

void Intervals::unionAll() {
	std::sort(_points.begin(), _points.end(), [](pair<bool, double> a, pair<bool, double> b) {
		if (a.second != b.second)
			return a.second < b.second;
		if (a.first != b.first)
			return a.first;
		return false;
	});

	size_t w = 0;
	size_t depth = 0;
	for (size_t r = 0; r < _points.size(); r++) {
		if (depth == 0 && _points[r].first)
			_points[w++] = _points[r];
		depth += _points[r].first ? 1 : -1;
		if (depth == 0 && !_points[r].first)
			_points[w++] = _points[r];
	}
	_points.resize(w);
}

int Classify(const xmesh3& m, const segment3& s, const aabb4& box) {
	if (!Intersects(box, aabb4(s)))
		return +1;

	Intervals intervals;
	for (const face& f : m) {
		vector<pair<double, double>> intersections;
		int c = Classify(f, s, &intersections);
		if (c == -1)
			return -1;
		if (c == 0)
			for (pair<double, double> p : intersections)
				intervals.add(p.first, p.second);
	}

	if (intervals.size() == 0)
		return +1;

	intervals.unionAll();

	// check all mid points between travels if they are inside mesh
	for (size_t i = 1; i < intervals.size(); i++) {
		double t = (intervals[i - 1].second + intervals[i].first) / 2;
		int c = Classify(m, s.linear(t), box);
		assert(c != 0);
		if (c == -1)
			return -1;
	}
	return 0;
}

static int ClassifyFaceEdgePairs(const xmesh3& ma, const xmesh3& mb, const aabb4& boxb) {
	int result = 1;
	for (const face& f : ma) {
		for (const segment3& e : Edges(f)) {
			int c = Classify(mb, e, boxb);
			if (c == 0)
				result = 0;
			if (c == -1)
				return -1;
		}
	}
	return result;
}

// This will be a ground truth function. Slow and accurate.
// TODO identify all contacts
int Classify(const xmesh3& ma, const xmesh3& mb) {
	aabb4 boxa(ma.vertices()), boxb(mb.vertices());
	if (!Intersects(boxa, boxb))
		return +1;

	if (Classify(ma, mb.vertices()[0], boxa) == -1 || Classify(mb, ma.vertices()[0], boxb) == -1)
		return -1;

	int ca = ClassifyFaceEdgePairs(ma, mb, boxb);
	if (ca == -1)
		return -1;
	int cb = ClassifyFaceEdgePairs(mb, ma, boxa);
	return min(ca, cb);
}
