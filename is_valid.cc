#include "is_valid.h"
#include "union_find.h"
#include "properties.h"
#include "util.h"
#include "aabb.h"
#include "exception.h"
#include "ivec.h"
#include "primitives.h"

// from tesselate.cc
bool relate_abxo(ivec2 a, ivec2 b, ivec2 p, ivec2 q, long ab);

bool is_valid(const ipolygon2& poly) {
	auto n = poly.size();
	if (n < 3)
		return false;

	// all vertices must be unique
	for (auto i : range(n))
		for (auto j : range(i))
			if (poly[i] == poly[j])
				return false;

	// no self intersection
	for (auto i : range(n)) {
		ivec2 a = poly[i], b = poly[(i + 1) % n];
		long ab = edge_area(a, b);
		for (auto j : range(i)) {
			ivec2 p = poly[j], q = poly[(j + 1) % n];
			if (relate_abxo(a, b, p, q, ab))
				return false;
		}
	}
	return true;
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
bool edge_overlap(const itriangle3& p, const itriangle3& q) {
	if (!aabb(p).intersects(aabb(q)))
		return false;
	for (isegment3 pe : p.edges())
		for (isegment3 qe : q.edges())
			if (colinear(pe.a, qe.a, qe.b) && colinear(pe.b, qe.a, qe.b))
				return aabb(pe).intersects(aabb(qe));
	return false;
}

// ===============

bool contains_multiple_components(const imesh3& mesh) {
	std::vector<UnionFind> component(mesh.size());
	for (auto i : range(mesh.size()))
		for (auto j : range(i))
			if (edge_overlap(mesh[i], mesh[j]))
				component[i].merge(component[j]);
	for (auto& c : component)
		if (c.find() != component[0].find())
			return true;
	return false;
}

bool face_cuts_into_other_face(const imesh3& mesh) {
	for (const itriangle3& m : mesh) {
		auto box = aabb(m.a);
		lvec3 ma = m.a, mb = m.b, mc = m.c;
		lvec3 cross_cb = crossi(mc, mb);
		lvec3 cross_ba = crossi(mb, ma);
		lvec3 cross_ac = crossi(ma, mc);
		lvec3 sub_cb = subi(mc, mb);
		lvec3 sub_ac = subi(ma, mc);

		for (const itriangle3& n : mesh)
			if (aabb<ivec3>(n).intersects(box))
				for (isegment3 e : n.edges()) {
					lvec3 ea = e.a, eb = e.b;
					lvec3 d = subi(eb, ea);
					lvec3 n = crossi(d, eb);
					long s = doti(n, sub_cb);
					long t = doti(n, sub_ac);
					if (doti(d, cross_cb) > negi(s) && doti(d, cross_ac) > negi(t) && doti(d, cross_ba) > addi(s, t))
						return true;
				}
	}
	return false;
}

bool contains_open_edge(const imesh3& mesh) {
	std::unordered_set<isegment3> open_edges;
	for (auto i : range(mesh.size()))
		for (auto e : mesh[i].edges())
			if (open_edges.erase(e.reversed()) == 0)
				open_edges.insert(e);
	for (isegment3 e : open_edges)
		print("open edge: %s\n", e);
	return !open_edges.empty();
}

// assume all three points are colinear
// assume A and B are on the same side of P
static bool closer_to(ivec3 p, ivec3 a, ivec3 b) {
	assert(colinear(a, b, p));
	assert(aabb<ivec3>(p, a).intersects(aabb<ivec3>(p, b)));
	return (p.x <= a.x && a.x < b.x) || (p.y <= a.y && a.y < b.y) || (p.z <= a.z && a.z < b.z)
		|| (p.x >= a.x && a.x > b.x) || (p.y >= a.y && a.y > b.y) || (p.z >= a.z && a.z > b.z);
}

// how many times can you add d to a before it overflows
static uint add_count(int a, long d) {
	if (d == 0)
		return std::numeric_limits<uint>::max();
	if (d > 0)
		return uint(std::numeric_limits<int>::max() - a) / uint(d);
	return uint(a - std::numeric_limits<int>::min()) / uint(-d);
}

static void minimize(ivec3& a, lvec3 d) {
	uint x = add_count(a.x, -d.x);
	uint y = add_count(a.y, -d.y);
	uint z = add_count(a.z, -d.z);
	a -= d * (long)min({x, y, z});
}

static void maximize(ivec3& a, lvec3 d) {
	uint x = add_count(a.x, d.x);
	uint y = add_count(a.y, d.y);
	uint z = add_count(a.z, d.z);
	a += d * (long)min({x, y, z});
}

struct Point {
	ivec3 pos;
	ivec3 angle_pos;
	ivec3 angle_off;
	bool begin;
	bool ccw; // when looking in the direction of line, is positive side of triangle pointing CCW?
};

void format_e(std::string& s, std::string_view spec, const Point& p) {
	format_s(s, "{pos: %s, angle_pos: %s, angle_off: %s, begin: %s, ccw: %s}",
		p.pos, p.angle_pos, p.angle_off, p.begin, p.ccw);
}

// angle between planes [line,A] and [line,REF]
// returns angle in range [-PI, PI]
double compute_angle(isegment3 line, ivec3 pos, ivec3 off, dvec3 dir, dvec3 ref, dvec3 cross_dir_ref) {
	dvec3 da = subi(off, pos);
	dvec3 e = glm::normalize(da - dir * glm::dot(da, dir));
	double alpha = angle(e, ref);
	return glm::dot(e, cross_dir_ref) < 0 ? -alpha : alpha;
}

struct LineData {
	std::vector<Point> points;
	lvec3 dir;
};

bool is_sealed(const imesh3& mesh) {
	// Extract all intervals on all lines of all edges
	std::unordered_map<isegment3, LineData> lines;
	for (const itriangle3& f : mesh)
		for (int i : range(3)) {
			ivec3 a = f[i], b = f[(i + 1) % 3], c = f[(i + 2) % 3];
			lvec3 d = (lvec3)b - (lvec3)a;
			d /= gcd(d);
			bool flip = d.x < 0 || (d.x == 0 && d.y < 0) || (d.x == 0 && d.y == 0 && d.z < 0);
			if (flip)
				d = -d;
			isegment3 e(a, b);
			minimize(e.a, d);
			maximize(e.b, d);
			auto& data = lines[e];
			data.points.push_back({a, flip ? a : b, c, !flip, !flip});
			data.points.push_back({b, flip ? a : b, c, flip, !flip});
			data.dir = d;
		}
	for (auto& [line, data] : lines) {
		ivec3 line_a = line.a;
		std::sort(data.points.begin(), data.points.end(), [&line_a](const Point& a, const Point& b) {
			return closer_to(line_a, a.pos, b.pos) || (a.pos == b.pos && !a.begin && b.begin);
		});
	}

	// Verify all edges are sealed and no seal is penetrating any other seal
	for (const auto& [line, data] : lines) {
		// coordinate system to compute angles in!
		dvec3 d = glm::normalize((dvec3)data.dir);
		dvec3 r = glm::normalize(any_normal(d));
		dvec3 dr = glm::normalize(glm::cross(d, r));

		auto less = [](std::pair<double, bool> a, std::pair<double, bool> b) {
			return a.first < b.first;
		};
		std::multiset<std::pair<double, bool>, decltype(less)> angles(less);
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

bool are_coplanar_faces_overlapping(const itriangle3& p, const itriangle3& q, lvec3 normal) {
	for (auto i : range(3)) {
		ivec3 ma = p[i];
		ivec3 mb = p[(i + 1) % 3];
		ivec3 mc = addi(ma, normal);
		ivec3 c = p[(i + 2) % 3];

		lvec3 normal = normali(ma, mb, mc);
		long pc = doti(subi(p.c, ma), normal);
		bool outside = true;
		for (ivec3 v : q) {
			long qa = doti(subi(v, ma), normal);
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

bool contains_overlapping_faces(const imesh3& mesh) {
	// are there two faces whose intersection is 2d?
	for (auto i : range(mesh.size()))
		for (auto j : range(i)) {
			const itriangle3& p = mesh[i], q = mesh[j];
			lvec3 normal = normali(p.a, p.b, p.c);
			if (doti(subi(q.a, p.a), normal) != 0)
				continue;
			if (doti(subi(q.b, p.a), normal) != 0)
				continue;
			if (doti(subi(q.c, p.a), normal) != 0)
				continue;
			if (are_coplanar_faces_overlapping(p, q, normal))
				return false;
		}
	return false;
}

Validity is_valid(const imesh3& mesh) {
	if (mesh.size() < 4)
		return Validity::TooFewFaces;

	// All faces must be valid
	for (const itriangle3& f : mesh)
		if (colinear(f.a, f.b, f.c))
			return Validity::InvalidFace;

	if (contains_multiple_components(mesh))
		return Validity::SeparateComponents;

	if (face_cuts_into_other_face(mesh))
		return Validity::SelfIntersection;

	if (contains_overlapping_faces(mesh))
		return Validity::OverlappingFaces;

	if (!is_sealed(mesh))
		return Validity::NotSealed;

	// check if all faces are oriented outside
	if (signed_volume_mul6(mesh) < 0)
		return Validity::Inverted;

	// Strict stronger condition than is_sealed
	//if (contains_open_edge(mesh))
	//	return Validity::OpenEdge;

	return Validity::OK;
}

void make_valid(imesh3& mesh) {
	THROW(not_implemented);
}
