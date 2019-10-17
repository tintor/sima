#include <core/exception.h>
#include <core/range.h>
#include <core/interval.h>

#include <geom/convex_body.h>
#include <geom/convex_hull.h>
#include <geom/pose.h>
#include <geom/sphere.h>
#include <geom/segment.h>

inline double hmin(double4 a) {
	// TODO auto b = vmin(a, vshuffle(a, a, 2, 3, 0, 1));
	return min(min(a.x, a.y), min(a.z, a.w));
}

inline double hmax(double4 a) {
	// TODO auto b = vmin(a, vshuffle(a, a, 2, 3, 0, 1));
	return max(max(a.x, a.y), max(a.z, a.w));
}

// pads by repeating the first element
static void Convert(const vector<double3>& a, vector<xyz4>& b) {
	b.resize((a.size() + 3) / 4);
	auto pad = [&](size_t i) { return (i < a.size()) ? a[i] : a[0]; };
	for (size_t i = 0; i < b.size(); i++) {
		double3 a0 = a[i * 4];
		double3 a1 = pad(i * 4 + 1);
		double3 a2 = pad(i * 4 + 2);
		double3 a3 = pad(i * 4 + 3);
		b[i].x = {a0.x, a1.x, a2.x, a3.x};
		b[i].y = {a0.y, a1.y, a2.y, a3.y};
		b[i].z = {a0.z, a1.z, a2.z, a3.z};
	}
}

inline interval<double> ProjectToAxis(const vector<xyz4>& convex, double3 axis) {
	const xyz4& e0 = convex[0];
	double4 smin = e0.x * axis.x + e0.y * axis.y + e0.z * axis.z;
	double4 smax = smin;

	for (size_t i = 1; i < convex.size(); i++) {
		const xyz4& e = convex[i];
		double4 s = e.x * axis.x + e.y * axis.y + e.z * axis.z;
		smin = vmin(smin, s);
		smax = vmax(smax, s);
	}
	return {hmin(smin), hmax(smax)};
}

inline interval<double> ProjectToAxis(const vector<double3>& convex, double3 axis) {
	interval<double> result;
	for (double3 p : convex) {
		result.add(dot(p, axis));
	}
	return result;
}

// direction is ignored
inline bool ContainsSimilarAxis(const vector<double3>& v, double3 a) {
	for (double3 e : v) {
		double t = angle(e, a);
		if (t <= Tolerance || PI - t <= Tolerance)
			return true;
	}
	return false;
}

cmesh3 GenerateConvexMesh(const vector<double3>& vertices) {
	cmesh3 m;
	Convert(vertices, m.vertices4);
	convex_hull(vertices, m.mesh);
	m.vertices = vertices;

	// TODO sort faceAxis by face area in decreasing order
	unordered_map<segment3, double3, hash_t<segment3>> edges;
	for (const auto& f : m.mesh) {
		double3 a = compute_normal(f.a, f.b, f.c);
		if (!ContainsSimilarAxis(m.faceAxis, a))
			m.faceAxis.push_back(a);
		for (auto e : Edges(f))
			edges.insert({e, a});
	}

	// TODO sort edgeAxis by edge length in decreasing order
	for (const auto& p : edges) {
		segment3 e = p.first;
		// ignore opposite edge and edges between coplanar faces
		if (!lex_less(e.a, e.b) && angle(edges.at(segment3{e.b, e.a}), p.second) > Tolerance) {
			double3 a = normalize(e.a - e.b);
			if (!ContainsSimilarAxis(m.edgeAxis, a))
				m.edgeAxis.push_back(a);
		}
	}
	return m;
}

inline double3 AxisBetweenConvexMeshes(const cmesh3& ca, const cmesh3& cb) {
	double3 axisBest;
	double distBest = -std::numeric_limits<double>::infinity();

	auto check = [&](double3 axis) {
		auto ia = ProjectToAxis(ca.vertices4, axis);
		auto ib = ProjectToAxis(cb.vertices4, axis);
		double dist = max(ia.min - ib.max, ib.min - ia.max);
		if (dist > distBest) {
			distBest = dist;
			axisBest = axis;
		}
	};

	for (double3 axis : ca.faceAxis)
		check(axis);
	for (double3 axis : cb.faceAxis)
		check(axis);
	for (double3 axisA : cb.edgeAxis)
		for (double3 axisB : cb.edgeAxis) {
			double3 axis = cross(axisA, axisB);
			double axisLen = length(axis);
			if (axisLen >= Tolerance)
				check(axis / axisLen);
		}

	return axisBest;
}

inline void FilterVertices(const vector<double3>& vertices, double3 axis, double dist, vector<double3>& out) {
	out.clear();
	for (double3 v : vertices) {
		double d = dot(v, axis);
		if (dist - Tolerance <= d && d <= dist + Tolerance)
			out.push_back(v);
	}
}

inline void ConvexHullInPlace2(vector<double3>& points) {
    if (points.size() <= 1)
		return;

    size_t leastY = 0;
    for (size_t i = 1; i < points.size(); i++)
        if (points[i].y < points[leastY].y)
            leastY = i;
    swap(points[0], points[leastY]);

	auto ccw = [](double3 a, double3 b, double3 c) {
	    return -( (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) );
	};

    double3 pivot = points[0];
	sort(points.begin() + 1, points.end(), [&pivot, &ccw](double3 a, double3 b)  {
	    double order = ccw(pivot, a, b);
	    if (order == 0)
    	    return squared(pivot - a) < squared(pivot - b);
	    return order < 0;
	});

	// remove dups after sorting
	points.erase(
		std::unique(points.begin(), points.end(), static_cast<bool(*)(double3, double3)>(&equal)),
		points.end());
    if (points.size() <= 3)
		return;

	size_t hull = 3;
	//print("start hull %s, points %s\n", hull, points);
    for (size_t i = 3; i < points.size(); i++) {
        auto top = points[hull - 1];
        hull -= 1;
		//print("A hull %s, points %s\n", hull, points);
        while (ccw(points[hull - 1], top, points[i]) >= 0) {
            top = points[hull - 1];
            hull -= 1;
			//print("B hull %s, points %s\n", hull, points);
        }
		points[hull++] = top;
		//print("C hull %s, points %s\n", hull, points);
		points[hull++] = points[i];
		//print("D hull %s, points %s\n", hull, points);
    }
	//print("end hull %s, points %s\n", hull, points);
	points.resize(hull);
}

#include <geom/classify.h>
inline double signed_double_edge_area2(double3 a, double3 b) {
	return (a.x + b.x) * (a.y - b.y);
}

inline double signed_double_area2(double3 a, double3 b, double3 c) {
	return signed_double_edge_area2(a, b) + signed_double_edge_area2(b, c) + signed_double_edge_area2(c, a);
}

inline int ClassifySegmentRayHelper2(double3 ea, double3 eb, double3 p) {
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

	return (signed_double_area2(ea, eb, p) < 0) ? -1 : 1;
}

inline int Classify2(segment3 s, double3 v) {
	return Classify(segment2(double2{s.a.x, s.a.y}, double2{s.b.x, s.b.y}), double2{v.x, v.y});
}

inline bool PointInPolygon2(double3 p, cspan<double3> ring) {
	for (auto e : Edges(ring))
		if (Classify2(e, p) == 0)
			return true;
	if (ring.size() < 3)
		return false;

	bool result = false;
	double entrance = 0;
	size_t m = IndexOfMax(ring, [p](double3 v) { return abs(v.y - p.y); });
	for (size_t i = 0; i < ring.size(); i++) {
		double3 a = ring[(i + m) % ring.size()];
		double3 b = ring[(i + m + 1) % ring.size()];

		int c = ClassifySegmentRayHelper2(a, b, p);
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

inline bool StrictIntersect2(double3 pa, double3 pb, double3 qa, double3 qb, double3& out) {
	segment2 p(double2{pa.x, pa.y}, double2{pb.x, pb.y});
	segment2 q(double2{qa.x, qa.y}, double2{qb.x, qb.y});
	double2 t;
	// TODO specialize relate for X case
	if ('X' == relate(p, q, &t, nullptr)) {
		auto e = p.linear(t.x);
		out = {e.x, e.y, 0};
		return true;
	}
	return false;
}

// TODO there is a faster way!
static void ConvexIntersect2(vector<double3>& a, vector<double3>& b, vector<double3>& c) {
	c.clear();
	// TODO use faster PointInConvexPolygon instead!
	double3 v;
	for (segment3 ea : Edges(a))
		for (segment3 eb : Edges(b))
			if (StrictIntersect2(ea.a, ea.b, eb.a, eb.b, v))
				c.push_back(v);

/*	if (c.size() == 0) {
		for (double3 v : a)
			if (PointInPolygon2(v, b)) {
				c = a;
				return;
			}
		for (double3 v : b)
			if (PointInPolygon2(v, a)) {
				c = b;
				return;
			}
		return;
	}*/

	for (double3 v : a)
		if (PointInPolygon2(v, b))
			c.push_back(v);
	for (double3 v : b)
		if (PointInPolygon2(v, a))
			c.push_back(v);

	// TODO this is just for ordering points, so a little cheaper algorithms could be used maybe?
	ConvexHullInPlace2(c);
}

static bool ProcessContacts(
	interval<double> ia,
	interval<double> ib,
	const cmesh3& ca,
	const cmesh3& cb,
	double3& normal,
	vector<double3>& contacts,
	vector<double3>& workA,
	vector<double3>& workB) {
	if ((ia.min + ia.max) > (ib.min + ib.max)) {
		normal = -normal;
		swap(ia.min, ia.max);
		swap(ib.min, ib.max);
		ia.min = -ia.min;
		ia.max = -ia.max;
		ib.min = -ib.min;
		ib.max = -ib.max;
	}

	// TODO find smallest abs component of normal and project cheaply by dropping it! (there will be some distortion)

	// TODO fast path for v/f: maybe just precompute edge planes and check them against 3d vertex without projection
	// TODO fast path for f/f: if all vertices of smaller face are inside edge planes of larger face just return smaller face
	// TODO fast path for f/e: if all vertices of edge are inside edge planes just return the edge
	// TODO fast path for e/e: run nearest(segment, segment) and return intersection (if only one point)

	assert(abs(ia.max - ib.min) <= Tolerance);
	double m = (ia.max + ib.min) / 2;
	// TODO if normal comes from one of the faces then pull out its vertices without filtering
	// TODO if normal comes from one of the edges then pull out its vertices without filtering
	FilterVertices(ca.vertices, normal, m, workA);
	FilterVertices(cb.vertices, normal, m, workB);

	double3 j = normalize(any_normal(normal));
	double3 k = cross(j, normal);
	for (double3& v : workA)
		v = {dot(v, j), dot(v, k), 0};
	for (double3& v : workB)
		v = {dot(v, j), dot(v, k), 0};

	// TODO instead of convex hull, check first three vertices and flip if clockwise
	ConvexHullInPlace2(workA);
	ConvexHullInPlace2(workB);

	ConvexIntersect2(workA, workB, contacts);
	if (contacts.size() == 0)
		return false;

	// TODO if vertex is original then reuse it!
	// TODO there should be a cheaper formula for unproject: origin + v.x * J + v.y * K
	// origin = solve_linear_row(j, k, normal, double3{0, 0, m}) ?
	for (double3& v : contacts)
		v = solve_linear_row(j, k, normal, double3{v.x, v.y, m});
	return true;
}

int ClassifyConvexConvex(
	const cmesh3& ca,
	const cmesh3& cb,
	bool hasNormal,
	double3& normal,
	vector<double3>& contacts,
	double& overlap,
	vector<double3>& workA,
	vector<double3>& workB) {

	int result = 0;
	overlap = std::numeric_limits<double>::infinity();
	double3 initialNormal = normal;

	auto check = [&](double3 axis) {
		if (hasNormal && equal(axis, initialNormal))
			return false;
		// TODO it is possible to pre-compute lower and upper bounds for an object for all axises to speed up disjoint check?
		// TODO maybe it helps to deduplicate axis list
		auto ia = ProjectToAxis(ca.vertices4, axis);
		auto ib = ProjectToAxis(cb.vertices4, axis);
		if (ia.max < ib.min - Tolerance || ib.max < ia.min - Tolerance) {
			normal = axis;
			result = 1;
			return true;
		}

		double axisOverlap = min(ia.max - ib.min, ib.max - ia.min);
		if (axisOverlap > Tolerance) {
			if (axisOverlap < overlap) {
				overlap = axisOverlap;
				normal = axis;
			}
			return false;
		}

		normal = axis;
		return ProcessContacts(ia, ib, ca, cb, normal, contacts, workA, workB);
	};

	if (hasNormal && check(normal))
		return result;

	// TODO object with the face with largest are goes first!
	for (double3 axis : ca.faceAxis)
		if (check(axis))
			return result;

	for (double3 axis : cb.faceAxis)
		if (check(axis))
			return result;

	// TODO better join iteration order (join longer edges before shorter edges)
	for (double3 axisA : cb.edgeAxis)
		for (double3 axisB : cb.edgeAxis) {
			double3 axis = cross(axisA, axisB);
			double axisLen = length(axis);
			if (axisLen >= Tolerance)
				if (check(axis / axisLen))
					return result;
		}

	return -1;
}
