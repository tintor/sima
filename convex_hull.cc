#include "convex_hull.h"
#include "is_valid.h"
#include "ivec.h"
#include "properties.h"
#include "aabb.h"
#include "util.h"
#include "primitives.h"

inline dvec3 dvec(ivec3 a) { return vconvert(a, dvec3); }

bool is_convex(const imesh3& mesh) {
	// TODO set might be unncessary complication
	// Extract all unique vertices
	std::unordered_set<ivec3, std::hash<ivec3>, equal_t<ivec3>> vertices;
	for (auto face : mesh)
		for (auto vertex : face)
			vertices.insert(vertex);
	// Check if every face has all vertices on its negative side
	for (auto face : mesh) {
		dvec3 n = normal(dvec(face.a), dvec(face.b), dvec(face.c));
		double d = dot(n, dvec(face.a));
		for (auto v : vertices)
			if (dot(n, dvec(v)) > d)
				return false;
	}
	return true;
}

	// TODO measure current convex hull for (rotated) bunny

	// SIMD scan
	// split vertices in 3 register streams
	// split 8 params in 3 registers
	// each param register does one 256bit multiply then add 3 registers together (and find min and max value)

	// 256 vectors
	// 8 directions (24 params)
	// 256bit = 8 floats
	// 3x 256 registers for 8 directions (register for each of X, Y, Z)

	// ideas:
	// - SIMD scan to find AABB
	// - find 4 unique and non-complanar points amoth those that are on AABB
	// - keep adding other AABB points to the hull (as they are definitely on it)
	// - (not needed after sort) fast SIMD scans to find extreme points (ie. points on the final hull) for 8 random directions
	// - construct initial hull from the collected points
	// - keep planes for hull faces (face can be more than triangle)
	// - SIMD scan through all points (discard points already in hull, assign remaining points to single face  that they are in front of face with max signed distance) some vertices will be discarded, others will be paritioned
	// - find the furthest vertex among all and add it next -> update distances of vertices that were in

// TODO maybe precompute this in separate array
uint boundary_dist(aabb<ivec3> box, ivec3 v) {
	std::array<uint, 3> m;
	for (auto i : range(3)) {
		int a = box.min[i], b = box.max[i];
		m[i] = (v[i] >= (a + b) / 2) ? (b - v[i]) : (v[i] - a);
	}
	return min(m[0], m[1], m[2]);
}

// convex_hull algorithms:
// A) naive:
//    - O(n) find xmin and xmax point for A and B
//    - O(n) find max distance from line AB for C
//    - O(n) find max distance from plane ABC
//    - O(n*h) add remaining points one by one to hull
//
// B) quick hull:
//    - find min and max points on each of main 3 axes (furthest ones are A and B)
//    - find max distance from line AB for C
//    - find max distance from plane ABC
//    - partition remaining points according to which plane they belong (or inside)
//    - process each partition one by one (furthest point first)
//    - after hull is updated, process all points in front of updated faces?
//
// C) scan these 13 "primary" axes using SIMD
// +00
// 0+0
// 00+
// 0++
// +0+
// ++0
// +++
// -++
// +-+
// ++-
// +-0
// +0-
// 0+-
// collect all vertices on each of the 26 "sides" and potentialy compute a 2d convex hull for them
// all of those 2d convex hulls will be polygonal faces on the big 3d hull
// How to combine them in initial hull?
//
// D) keep rounding point coordinates until we are left with small number of unique points in 3d grid (at least 4)
//    build hull from those initial 4, and then keep unrounding points slowly and expanding hull
//    idea is to build hull of simplified point set first, which is a "good" approximation of the point set
//    (problem, one dimension could be thinner than other two, so after a few rounding we end up with 2d case)

void convex_hull_fast(array_cptr<ivec3> points, imesh3& hull) {
	hull.clear();
	if (points.size() < 4)
		return;

	// 1) compute aabb
	aabb<ivec3> box(points);

	// 2) sort points by reverse manhattan distance from the AABB (after manhattan dist, use lex order)
	std::unique_ptr<uint[]> order(new uint[points.size()]);
	for (uint i = 0; i < points.size(); i++)
		order[i] = i;
	const ivec3* pts = points.begin();
	std::sort(order.get(), order.get() + points.size(), [pts, &box] (uint a, uint b) {
		ivec3 A = pts[a], B = pts[b];
		uint da = boundary_dist(box, A);
	    uint db = boundary_dist(box, B);
		if (da < db) return true;
		if (da > db) return false;
		if (A.x < B.x) return true;
		if (A.x > B.x) return false;
		if (A.y < B.y) return true;
		if (A.y > B.y) return false;
		return A.z < B.z;
	});

	// 3) find initial 4 hull points
	// either: first four points in order will be non-complanar -> build hull
	//      OR some points will be duplicates
	//      OR some points will be colinear / coplanar with others
	//   try to maximize

	// 4a) add points one by one in increasing distance from hull order

	// 4b) quick-hull: partition remaining points into
	THROW(not_implemented);
}

void convex_hull(array_cptr<ivec3> points, imesh3& hull) {
	hull.clear();
	if (points.size() < 4)
		return;

	// First two points (A and B) on the hull (extremes on X axis)
	ivec3 a = points[0], b = points[0];
	for (auto p : points) {
		if (p.x < a.x)
			a = p;
		if (p.x > b.x)
			b = p;
	}
	if (equal(a, b))
		return;

	// Third point C on the hull (furthest from line AB)
	int128 max_dist2 = 0;
	ivec3 c;
	for (auto p : points) {
		lvec3 cross = crossi(subi(p, a), subi(p, b));
		int128 dist2 = squaredi(vconvert(cross, llvec3));
		if (dist2 > max_dist2) {
			c = p;
			max_dist2 = dist2;
		}
	}
	if (max_dist2 == 0)
		return;

	// Fourth point D on the hull (furthest from plane ABC)
	int128 max_dist = 0;
	ivec3 d;
	lvec3 normal = normali(a, b, c);
	long dd = doti(normal, vconvert(a, lvec3));
	for (auto p : points) {
		int128 dist = subi(doti(vconvert(normal, llvec3), vconvert(p, llvec3)), dd);
		if (std::abs(dist) > std::abs(max_dist)) {
			d = p;
			max_dist = dist;
		}
	}
	if (max_dist == 0)
		return;

	// Construct initial tetrahedron hull,
	// All faces are oriented facing outside with right hand rule.
	hull.reserve(4);
	if (max_dist < 0) {
		hull.emplace_back(a, b, c);
		hull.emplace_back(b, a, d);
		hull.emplace_back(c, b, d);
		hull.emplace_back(a, c, d);
	} else {
		hull.emplace_back(c, b, a);
		hull.emplace_back(a, b, d);
		hull.emplace_back(b, c, d);
		hull.emplace_back(c, a, d);
	}

	// Expand hull to include all remaining points outside of it
	std::unordered_set<isegment3> open_edges;
	uint i = 0;
	for (auto p : points) {
		if (points.size() >= 50000 && i++ % (1 << 14) == 0)
			print("%.2f hull_vertices=%s hull_volume=%s\n", double(i) / points.size(), hull.size(), volume(hull));
		if (equal(p, a) || equal(p, b) || equal(p, c) || equal(p, d))
			continue;
		// Remove faces on hull covered by new vertex
		open_edges.clear();
		for (size_t i = 0; i < hull.size(); i++) {
			const itriangle3& f = hull[i];
			// Skip if P is not in front of face F
			lvec3 n = normali(f.a, f.b, f.c);
			if (doti(vconvert(n, llvec3), vconvert(subi(p, f.a), llvec3)) <= 0)
				continue;
			// Add edges of removed face to open_edges
			for (auto e : f.edges()) {
				// If two faces that share an edge are both removed,
				// then their shared edge isn't open anymore.
				if (open_edges.erase(e.reversed()) == 0)
					open_edges.insert(e);
			}
			// Remove face
			hull[i] = hull[hull.size() - 1];
			hull.resize(hull.size() - 1);
			i -= 1;
		}
		// For each open edge create a face that connects it with P
		for (auto e : open_edges)
			hull.emplace_back(e.a, e.b, p);
	}

	assert(is_valid(hull) == Validity::OK);
	assert(is_convex(hull));
}
