#include "classify.h"
#include "exception.h"
#include "edges.h"
#include "project.h"
#include "properties.h"

int Classify(const xpolygon2& f, ray2 s) {
	THROW(not_implemented);
}

static int ClassifySegmentRayHelper(double2 ea, double2 eb, double2 p) {
	if (ea.y > eb.y)
		swap(ea, eb);

	int sa = Sign(ea.y - p.y);
	int sb = Sign(eb.y - p.y);
	if (sa > 0 || sb < 0)
		return 1;

	if (sa == 0 && sb == 0)
		return (p.x > max(ea.x, eb.x) + Tolerance) ? 1 : 0;
	if (sa == 0)
		return (p.x > ea.x + Tolerance) ? 1 : 0;
	if (sb == 0)
		return (p.x > eb.x + Tolerance) ? 1 : 0;

	return (signed_double_area(ea, eb, p) < 0) ? -1 : 1;
}

static bool PointInPolygon(double2 p, span<const double2> ring) {
	bool result = false;
	double entrance = 0;
	size_t m = IndexOfMax(ring, [p](double2 v) { return abs(v.y - p.y); });
	for (size_t i = 0; i < ring.size(); i++) {
		double2 a = ring[(i + m) % ring.size()];
		double2 b = ring[(i + m + 1) % ring.size()];

		int c = ClassifySegmentRayHelper(a, b, p);
		if (c == -1)
			result ^= 1;
		if (c == 0) {
			bool za = abs(a.y - p.y) <= Tolerance;
			bool zb = abs(b.y - p.y) <= Tolerance;
			if (zb && !za)
				entrance = a.y - p.y;
			if (za && !zb) {
				if (entrance * (b.y - p.y) < 0)
					result ^= 1;
				entrance = 0;
			}
		}
	}
	return result;
}

template<typename Polygon2>
int TClassify(const Polygon2& a, double2 p, aabb2 box) {
	if (!box.intersects(p))
		return +1;

	for (auto e : Edges(a))
		if (Classify(e, p) == 0)
			return 0;

	int result = 1;
	for (auto ring : Rings(a))
		if (PointInPolygon(p, ring))
			result = -result;
	return result;
}

int Classify(const polygon2& a, double2 p, aabb2 box) { return TClassify(a, p, box); }
int Classify(const xpolygon2& a, double2 p, aabb2 box) { return TClassify(a, p, box); }

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
int TClassify(const Polygon2& f, segment2 s, aabb2 box, vector<pair<double, double>>* intersections = nullptr) {
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

int Classify(const polygon2& f, segment2 s, aabb2 box, vector<dpair>* intersections) {
	return TClassify(f, s, box, intersections);
}
int Classify(const xpolygon2& f, segment2 s, aabb2 box, vector<dpair>* intersections) {
	return TClassify(f, s, box, intersections);
}

template<typename Polygon2>
int TClassify(const Polygon2& a, const Polygon2& b) {
	aabb2 va = Box(a), vb = Box(b);
	if (!Intersects(va, vb))
		return +1;

	if (Contains(va, vb) && Classify(a, AnyVertex(b), va) < 0)
		return -1;
	if (Contains(vb, va) && Classify(b, AnyVertex(a), vb) < 0)
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

int Classify(const xpolygon2& a, const xpolygon2& b) { return TClassify(a, b); }
int Classify(const polygon2& a, const polygon2& b) { return TClassify(a, b); }

int Classify(const face& f, double4 v, const aabb4& box) {
	if (!Intersects(box, aabb4(v)))
		return +1;

	plane p = f.plane();
	double d = p.distance(v);
	if (abs(d) > Tolerance)
		return +1;

	int axis = ProjectionAxis(f);
	int c = Classify(Project(f, axis), Project(v - p.normal() * d, axis));
	return (c == 1) ? 1 : 0;
}

int Classify(const face& f, const ray3& s, const aabb4& box, double* travel) {
	if (!Intersects(box, aabb4(s.origin, s.infinity())))
		return +1;

	plane p = f.plane();
	double d = p.distance(s.origin);
	double t = -d / p.distance(s.unit_dir);
	if (travel)
		*travel = t;
	if (std::isfinite(t)) {
		if (t < -Tolerance)
			return +1;
		return Classify(f, s.linear(t), box);
	}

	if (abs(d) > Tolerance)
		return +1;

	int axis = ProjectionAxis(f);
	int c = Classify(Project(f, axis), Project(s, axis));
	return (c == 1) ? 1 : 0;
}

pair<int, int> ClassifyDoubleSided(const xpolygon2& poly, const ray2& ray) {
	THROW(not_implemented);
}

// Computes Classify for both s and -s at once.
pair<int, int> ClassifyDoubleSided(const face& f, const ray3& s, const aabb4& box) {
	if (!Intersects(box, aabb4(s.origin, s.infinity())) && !Intersects(box, aabb4(s.origin, -s.infinity())))
		return {1, 1};

	plane p = f.plane();
	double d = p.distance(s.origin);
	double t = -d / p.distance(s.unit_dir);
	if (std::isfinite(t)) {
		// ray intersects plane in point
		int c = Classify(f, s.linear(t), box);
		if (c == 1)
			return {1, 1};
		if (t > Tolerance)
			return {c, 1};
		if (t < -Tolerance)
			return {1, c};
		return {c, c};
	}

	// ray is parallel with plane
	if (d > Tolerance || d < -Tolerance)
		return {1, 1};

	// ray is coplanar with plane
	int axis = ProjectionAxis(f);
	auto c = ClassifyDoubleSided(Project(f, axis), Project(s, axis));
	if (c.first == -1)
		c.first = 0;
	if (c.second == -1)
		c.second = 0;
	return c;
}

int Classify(plane p, const segment3& s) {
	return Sign(p, s.a) * Sign(p, s.b);
}

// [intersections] represent relative travels along segment S
int Classify(const face& f, const segment3& s, bool reverse, vector<dpair>* intersections, vector<Contact>* contacts) {
	if (!Intersects(Box(f.vertices()), aabb4(s)))
		return +1;

	plane p = f.plane();
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
		int c = Classify(f2, Project(s, axis), Box(f2), intersections);
		if (contacts) {
			if (c == 0) {
				// edge V vertex/edge contact
				Contact contact;
				contact.type = Contact::Type::Point;
				auto i = intersections->back();
				contact.sa = s.linear(i.first);
				contact.normal = p.normal(); // TODO normal depends on two edges
				contacts->push_back(contact);
			}
			if (c < 0) {
				// edge V foce contact
				Contact contact;
				contact.type = Contact::Type::Segment;
				auto i = intersections->back();
				contact.sa = s.linear(i.first);
				contact.sb = s.linear(i.second);
				contact.normal = reverse ? -p.normal() : p.normal();
				contacts->push_back(contact);
			}
		}
		return c;
	}

	if (ca == 0 || cb == 0) {
		// segment touches the plane
		auto v2 = Project((ca == 0) ? s.a : s.b, axis);
		int c = Classify(f2, v2);
		if (c == 1)
			return 1;
		if (contacts) {
			if (c < 0) {
				// vertex V face contact
				Contact contact;
				contact.type = Contact::Type::Point;
				contact.sa = (ca == 0) ? s.a : s.b;
				contact.normal = reverse ? -p.normal() : p.normal();
				contacts->push_back(contact);
			}
			if (c == 0) {
				// vertex V edge contact
				Contact contact;
				contact.type = Contact::Type::Point;
				contact.sa = (ca == 0) ? s.a : s.b;
				contact.normal = p.normal(); // TODO normal?
				contacts->push_back(contact);
			}
		}
		if (intersections) {
			double t = (ca == 0) ? 0 : 1;
			intersections->emplace_back(t, t);
		}
		return 0;
	}

	// segment intersects the plane
	double t = da / (da - db);
	auto v = s.linear(t);
	int c = Classify(f2, Project(v, axis));
	if (c == 1)
		return 1;
	if (c == 0 && contacts) {
		// edge V edge contact
		Contact contact;
		contact.type = Contact::Type::Point;
		contact.sa = s.linear(t);
		//TODO contact.normal = cross(s.b - s.a, ?);
		contacts->push_back(contact);
	}
	if (intersections)
		intersections->emplace_back(t, t);
	return c;
}

int Classify(const xmesh3& m, double4 p, const aabb4& box) {
	if (!box.intersects(p))
		return +1;

	// Check if the point is on the boundary
	for (const face& f : m)
		if (Classify(f, p, box) == 0)
			return 0;

	// Shoot a ray (in both directions at once) and count crossings
	// (if both rays hits any edge or vertex, shoot a new ray and repeat)
	std::default_random_engine rnd;
	std::uniform_int_distribution dist(0, 2);
	rnd.seed(0);
	double4 dirs[3] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}};
	for (size_t i = 0; ; i += 1) {
		// ray direction vector will have 1 or 2 components zero to make its bounding box smaller
		double4 d;
		if (i < 3) {
			d = dirs[i];
		} else {
			while (true) {
				d = uniform_dir3(rnd);
				d[dist(rnd)] = 0;
				if (squared(d) > squared(0.01))
					break;
			}
		}
		ray3 ray(p, p + d);
		pair<int, int> result = {1, 1};
		for (const face& f : m) {
			auto c = ClassifyDoubleSided(f, ray, box);
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

int Classify(const xmesh3& m, const segment3& s, const aabb4& box, bool reverse, vector<Contact>* contacts) {
	if (!Intersects(box, aabb4(s)))
		return +1;

	Intervals intervals;
	for (const face& f : m) {
		vector<dpair> intersections;
		int c = Classify(f, s, reverse, &intersections, contacts);
		if (c == -1)
			return -1;
		if (c == 0)
			for (dpair p : intersections)
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

static int ClassifyFaceEdgePairs(const xmesh3& ma, const xmesh3& mb, const aabb4& boxb, bool reverse, vector<Contact>* contacts) {
	int result = 1;
	for (const face& f : ma)
		for (const segment3& e : Edges(f)) {
			// Since meshes are assumed to be closed, each edge will appear twice (in opposited directions)
			if (e.a.x > e.b.x)
				continue;

			int c = Classify(mb, e, boxb, reverse, contacts);
			if (c == 0)
				result = 0;
			if (c == -1)
				return -1;
		}
	return result;
}

bool InteriorIntersects(segment3 p, segment3 q, double4* common) {
	double4 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	double aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr double tiny = 1e-8;

	double d = aa * bb - ab * ab;
	double s = ab * bc - bb * ac;
	double t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0)) {
		double4 ma = p.a + A * (s / d);
		double4 mb = q.a + B * (t / d);
		if (Equals(ma, mb) && !Equals(ma, p.a) && !Equals(ma, p.b) && !Equals(mb, q.a) && !Equals(mb, q.b)) {
			*common = (ma + mb) / 2;
			return true;
		}
	}
	return false;
}

// This will be a ground truth function. Slow and accurate.
// TODO identify all contacts
int Classify(const xmesh3& ma, const xmesh3& mb, vector<Contact>* contacts) {
	auto va = Box(ma.vertices());
	auto vb = Box(mb.vertices());
	if (!Intersects(va, vb))
		return +1;

	if (Contains(va, vb) && Classify(ma, mb.vertices()[0], va) < 0)
		return -1;
	if (Contains(vb, va) && Classify(mb, ma.vertices()[0], vb) < 0)
		return -1;

	int ca = ClassifyFaceEdgePairs(ma, mb, vb, false, contacts);
	if (ca == -1)
		return -1;
	int cb = ClassifyFaceEdgePairs(mb, ma, va, true, contacts);
	if (cb == -1)
		return -1;
	int c = min(ca, cb);

	if (contacts && c == 0) {
		// edge cross edge contacts
		double4 common;
		for (auto [ea, na] : ma.uniqueEdges())
			for (auto [eb, nb] : mb.uniqueEdges())
				if (InteriorIntersects(ea, eb, &common)) {
					Contact contact;
					contact.type = Contact::Type::Point;
					contact.sa = common;
					contact.normal = normalize(cross(ea.b - ea.a, eb.b - eb.a)); // TODO direction
					contacts->push_back(contact);
				}

		// face touch vertex contacts
		// (includes vertex V vertex/edge cases)
		for (const auto& f : ma)
			for (auto v : mb.uniqueVertices()) {
				double d = f.plane().distance(v);
				if (abs(d) <= Tolerance) {
					// TODO point in polygon check
				}
			}
		for (const auto& f : mb)
			for (auto v : ma.uniqueVertices()) {
				double d = f.plane().distance(v);
				if (abs(d) <= Tolerance) {
					// TODO point in polygon check
				}
			}
	}
	return c;
}
