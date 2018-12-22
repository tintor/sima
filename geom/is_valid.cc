#include <geom/is_valid.h>
#include <core/union_find.h>
#include "properties.h"
#include <core/util.h>
#include <geom/aabb.h>
#include <core/exception.h>
#include "primitives.h"
#include "project.h"
#include <geom/classify.h>

bool IsValidPolygon(cspan<double2> poly) {
	auto n = poly.size();
	if (n < 3)
		return false;

	// all vertices must be unique
	for (auto i : range(n))
		for (auto j : range(i))
			if (equal(poly[i], poly[j]))
				return false;

	// no self intersection
	for (auto i : range(n)) {
		double2 a = poly[i], b = poly[(i + 1) % n];
		double inv_ab_len = 1 / length(a - b);
		for (auto j : range(i)) {
			double2 p = poly[j], q = poly[(j + 1) % n];
			if (relate_abxo(segment2(a, b), segment2(p, q), inv_ab_len))
				return false;
		}
	}
	return true;
}

bool IsValid(const polygon2& poly) {
	return IsValidPolygon(poly);
}

bool IsValid(const xpolygon2& poly) {
	THROW(not_implemented);
}

// is triangle Q intersecting with plane defined by plane of triangle P
/*inline bool intersects_plane(const itriangle3& q, const itriangle3& plane) {
	// TODO is this overflow safe?
	auto n = cross(plane.b - plane.a, plane.c - plane.a);
	auto e = dot(n, plane.a);

	auto a = dot(n, q.a) - e;
	auto b = dot(n, q.b) - e;
	auto c = dot(n, q.a) - e;
	return (a <= 0 || b <= 0 || c <= 0) && (a >= 0 || b >= 0 || c >= 0);
}*/

// if vertex/edge touch only return false
bool edge_overlap(triangle3 p, triangle3 q) {
	if (!Intersects(aabb4(p), aabb4(q)))
		return false;
	for (segment3 pe : Edges(p))
		for (segment3 qe : Edges(q))
			if (colinear(pe.a, qe.a, qe.b) && colinear(pe.b, qe.a, qe.b))
				return Intersects(aabb4(pe), aabb4(qe));
	return false;
}

// ===============

bool contains_multiple_components(const mesh3& mesh) {
	vector<UnionFind> component(mesh.size());
	for (auto i : range(mesh.size()))
		for (auto j : range(i))
			if (edge_overlap(mesh[i], mesh[j]))
				component[i].merge(component[j]);
	for (auto& c : component)
		if (c.find() != component[0].find())
			return true;
	return false;
}

bool face_cuts_into_other_face(const mesh3& mesh) {
	for (const triangle3& m : mesh) {
		auto box = aabb4(m.a);
		double4 cross_cb = cross(m.c, m.b);
		double4 cross_ba = cross(m.b, m.a);
		double4 cross_ac = cross(m.a, m.c);
		double4 sub_cb = m.c - m.b;
		double4 sub_ac = m.a - m.c;

		for (triangle3 n : mesh)
			if (Intersects(aabb4(n), box))
				for (segment3 e : Edges(n)) {
					double4 d = e.b - e.a;
					double4 n = cross(d, e.b);
					long s = dot(n, sub_cb);
					long t = dot(n, sub_ac);
					if (dot(d, cross_cb) > -s && dot(d, cross_ac) > -t && dot(d, cross_ba) > s + t)
						return true;
				}
	}
	return false;
}

bool contains_open_edge(const mesh3& mesh) {
	unordered_set<segment3, hash_t<segment3>> open_edges;
	for (auto i : range(mesh.size()))
		for (auto e : mesh[i].edges())
			if (open_edges.erase(e.reversed()) == 0)
				open_edges.insert(e);
	for (segment3 e : open_edges)
		print("open edge: %s\n", e);
	return !open_edges.empty();
}

// assume all three points are colinear
// assume A and B are on the same side of P
static bool closer_to(double4 p, double4 a, double4 b) {
	assert(colinear(a, b, p));
	assert(aabb4(p, a).intersects(aabb4(p, b)));
	return (p.x <= a.x && a.x < b.x) || (p.y <= a.y && a.y < b.y) || (p.z <= a.z && a.z < b.z)
		|| (p.x >= a.x && a.x > b.x) || (p.y >= a.y && a.y > b.y) || (p.z >= a.z && a.z > b.z);
}

struct Point {
	double4 pos;
	double4 angle_pos;
	double4 angle_off;
	bool begin;
	bool ccw; // when looking in the direction of line, is positive side of triangle pointing CCW?
};

void format_e(string& s, string_view spec, const Point& p) {
	format_s(s, "{pos: %s, angle_pos: %s, angle_off: %s, begin: %s, ccw: %s}",
		p.pos, p.angle_pos, p.angle_off, p.begin, p.ccw);
}

// angle between planes [line,A] and [line,REF]
// returns angle in range [-PI, PI]
// TODO [line] is not used?
double compute_angle(line3 line, double4 pos, double4 off, double4 dir, double4 ref, double4 cross_dir_ref) {
	double4 da = off - pos;
	double4 e = normalize(da - dir * dot(da, dir));
	double alpha = angle(e, ref);
	return dot(e, cross_dir_ref) < 0 ? -alpha : alpha;
}

struct LineData {
	vector<Point> points;
	double4 dir;
};

bool IsSealed(const mesh3& mesh) {
	// Extract all intervals on all lines of all edges
	unordered_map<line3, LineData, hash_t<line3>> lines;
	for (const triangle3& f : mesh)
		for (int i : range(3)) {
			double4 a = f[i], b = f[(i + 1) % 3], c = f[(i + 2) % 3];
			double4 d = normalize(b - a);
			bool flip = d.x < 0 || (d.x == 0 && d.y < 0) || (d.x == 0 && d.y == 0 && d.z < 0);
			if (flip)
				d = -d;
			line3 e{a, b};
			auto& data = lines[e];
			data.points.push_back({a, flip ? a : b, c, !flip, !flip});
			data.points.push_back({b, flip ? a : b, c, flip, !flip});
			data.dir = d;
		}
	for (auto& [line, data] : lines) {
		double4 line_a = line.origin;
		std::sort(data.points.begin(), data.points.end(), [&line_a](const Point& a, const Point& b) {
			// TODO verify ordering logic!
			return closer_to(line_a, a.pos, b.pos) || (equal(a.pos, b.pos) && !a.begin && b.begin);
		});
	}

	// Verify all edges are sealed and no seal is penetrating any other seal
	for (const auto& [line, data] : lines) {
		// coordinate system to compute angles in!
		double4 d = data.dir;
		double4 r = normalize(any_normal(d));
		double4 dr = normalize(cross(d, r));

		auto less = [](pair<double, bool> a, pair<double, bool> b) {
			return a.first < b.first;
		};
		multiset<pair<double, bool>, decltype(less)> angles(less);
		for (auto a = data.points.begin(); a != data.points.end(); a++) {
			double angle = compute_angle(line, a->angle_pos, a->angle_off, d, r, dr);
			if (a->begin)
				angles.insert({angle, a->ccw});
			else
				angles.erase({angle, a->ccw});
			auto b = a + 1;
			if (a->begin && b != data.points.end() && !b->begin)
				for (auto p = angles.begin(); p != angles.end(); p++) {
					auto q = p;
					++q;
					if (q == angles.end())
						q = angles.begin();
					if (p->first == q->first || p->second == q->second)
						return false;
				}
		}
		assert(angles.empty());
	}
	return true;
}

bool are_coplanar_faces_overlapping(triangle3 p, triangle3 q, double4 normal) {
	for (auto i : range(3)) {
		double4 ma = p[i];
		double4 mb = p[(i + 1) % 3];
		double4 mc = ma + normal;
		double4 c = p[(i + 2) % 3];

		double4 normal = compute_normal(ma, mb, mc);
		double pc = dot(p.c - ma, normal); // TODO possible typo bug
		bool outside = true;
		for (double4 v : q) {
			double qa = dot(v - ma, normal);
			if ((pc > 0 && qa > 0) || (pc < 0 && qa < 0)) {
				outside = false;
				break;
			}
		}
		if (outside)
			return false;
	}
	return true;
}

bool contains_overlapping_faces(const mesh3& mesh) {
	// are there two faces whose intersection is 2d?
	for (auto i : range(mesh.size())) {
		auto p = mesh[i];
		double4 n = compute_normal(p.a, p.b, p.c);
		for (auto j : range(i)) {
			auto q = mesh[j];
			if (dot(q.a - p.a, n) != 0)
				continue;
			if (dot(q.b - p.a, n) != 0)
				continue;
			if (dot(q.c - p.a, n) != 0)
				continue;
			if (are_coplanar_faces_overlapping(p, q, n))
				return false;
		}
	}
	return false;
}

Validity IsValid(const mesh3& mesh) {
	if (mesh.size() < 4)
		return Validity::TooFewFaces;

	// All faces must be valid
	for (auto f : mesh)
		if (colinear(f[0], f[1], f[2]))
			return Validity::InvalidFace;

	if (contains_multiple_components(mesh))
		return Validity::SeparateComponents;

	if (face_cuts_into_other_face(mesh))
		return Validity::SelfIntersection;

	if (contains_overlapping_faces(mesh))
		return Validity::OverlappingFaces;

	if (!IsSealed(mesh))
		return Validity::NotSealed;

	// check if all faces are oriented outside
	if (SignedVolume(mesh) < 0)
		return Validity::Inverted;

	// Strict stronger condition than is_sealed
	//if (contains_open_edge(mesh))
	//	return Validity::OpenEdge;

	return Validity::OK;
}

bool AreCoplanarFacesOverlapping(face p, face q) {
	int axis = ProjectionAxis(p);
	xpolygon2 p2 = Project(p, axis);
	xpolygon2 q2 = Project(q, axis);
	// TODO issue this projection will stretch some axis more than others and skew tolerances!
	return Classify(p2, q2) == -1;
}

// +2 - disjoint positive
// +1 - touching from positive
//  0 - coplanar
// -1 - touching from negative
// -2 - disjoint negative
// -100 - crossing
int Classify(plane p, face f) {
	int count[3] = {0, 0, 0};
	for (point3 v : f.vertices())
		count[Sign(p, v) + 1] += 1;
	int neg = count[0];
	int zero = count[1];
	int pos = count[2];

	if (neg == 0 && pos > 0)
		return (zero > 0) ? +1 : +2;
	if (neg > 0 && pos == 0)
		return (zero > 0) ? -1 : -2;
	if (neg == 0 && zero == 0 && pos == 0)
		return 0;
	return -100;
}

bool AreCoplanar(plane p, face f) {
	for (point3 v : f.vertices())
		if (Sign(p, v) != 0)
			return false;
	return true;
}

bool ContainsOverlappingFaces(const xmesh3& mesh) {
	// are there two faces whose intersection is 2d?
	for (auto i : range(mesh.size())) {
		face p = mesh[i];
		plane pp = best_fit_plane(p.vertices());
		for (auto j : range(i)) {
			face q = mesh[j];
			if (AreCoplanar(pp, q) && AreCoplanarFacesOverlapping(p, q))
				return true;
		}
	}
	return false;
}

// each edge must have its reversed edge
bool IsSealed(const xmesh3& mesh) {
	unordered_multiset<segment3, hash_t<segment3>> edges;
	for (face f : mesh)
		for (auto ring : f)
			for (auto edge : Edges(ring))
				edges.insert(edge);
	for (const segment3& e : edges)
		if (edges.count(e) != edges.count(e.reversed()))
			return false;
	return true;
}

Validity IsValid(const xmesh3& mesh) {
	if (mesh.size() < 4)
		return Validity::TooFewFaces;

	// All faces must be valid
	for (face f : mesh) {
		if (!AreCoplanar(best_fit_plane(f.vertices()), f))
			return Validity::NonPlanarFace;
		xpolygon2 poly = Project(f, ProjectionAxis(f));
		if (!IsValid(poly))
			return Validity::InvalidFace;
	}

	/*if (contains_multiple_components(mesh))
		return Validity::SeparateComponents;

	if (face_cuts_into_other_face(mesh))
		return Validity::SelfIntersection;*/

	if (ContainsOverlappingFaces(mesh))
		return Validity::OverlappingFaces;

	if (!IsSealed(mesh))
		return Validity::NotSealed;

	// check if all faces are oriented outside
	if (SignedVolume(mesh) < 0)
		return Validity::Inverted;

	return Validity::OK;
}

void MakeValid(mesh3& mesh) {
	THROW(not_implemented);
}

void MakeValid(xmesh3& mesh) {
	THROW(not_implemented);
}
