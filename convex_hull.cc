#include "convex_hull.h"
#include <unordered_set>

bool is_convex(const imesh3& mesh) {
	// Extract all unique vertices
	std::unordered_set<ivec3> vertices;
	for (auto face : mesh)
		for (auto vertex : face)
			vertices.insert(vertex);
	// Check if every face has all vertices on its negative side
	for (auto face : mesh) {
		auto normal = cross(face.b - face.a, face.c - face.a);
		auto d = dot(normal, face.a);
		for (auto v : vertices)
			if (dot(normal, v) > d)
				return false;
	}
	return true;
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
	if (a == b)
		return;

	// Third point C on the hull (furthest from line AB)
	int128 max_dist2 = 0;
	ivec3 c;
	for (auto p : points) {
		auto dist2 = squared(cross(p - a, p - b));
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
	lvec3 normal = cross(b - a, c - a);
	long dd = dot(normal, a);
	for (auto p : points) {
		int128 dist = dot(normal, p) - dd;
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
	assert(is_convex(hull));

	// Expand hull to include all remaining points outside of it
	std::unordered_set<isegment3> open_edges;
	for (auto p : points) {
		if (p == a || p == b || p == c || p == d)
			continue;
		print("%s\n", hull);
		// Remove faces on hull covered by new vertex
		open_edges.clear();
		for (auto i : range(hull.size())) {
			const itriangle3& f = hull[i];
			// Skip if P is not in front of face F
			if (dot(cross(f.b - f.a, f.c - f.a), p - f.a) <= 0)
				continue;
			// Add edges of removed face to open_edges
			for (auto e : f.edges())
				// If two faces that share an edge are both removed,
				// then their shared edge isn't open anymore.
				if (open_edges.erase(e.reversed()) == 0)
					open_edges.insert(e);
			// Remove face
			hull[i] = hull[hull.size() - 1];
			hull.resize(hull.size() - 1);
			i -= 1;
		}
		// For each open edge create a face that connects it with P
		for (auto e : open_edges)
			hull.emplace_back(e.a, e.b, p);
		print("%s\n", hull);
		assert(is_convex(hull));
	}
}

