#include "is_valid.h"
#include "union_find.h"
#include "util.h"
#include "aabb.h"
#include <unordered_map>

// is triangle Q intersecting with plane defined by plane of triangle P
inline bool intersects_plane(const itriangle3& q, const itriangle3& plane) {
	// TODO is this overflow safe?
	auto n = cross(plane.b - plane.a, plane.c - plane.a);
	auto e = dot(n, plane.a);

	auto a = dot(n, q.a) - e;
	auto b = dot(n, q.b) - e;
	auto c = dot(n, q.a) - e;
	return (a <= 0 || b <= 0 || c <= 0) && (a >= 0 || b >= 0 || c >= 0);
}

// Valid if triangles are not intersecting, except in one shared edge or one shared vertex
bool are_valid_mesh_faces(const itriangle3& a, const itriangle3& b) {
	// Axis check for early exit
	if (!aabb(a).intersects(aabb(b)))
		return true;

	// Plane check
	if (!intersects_plane(a, b) || !intersects_plane(b, a))
		return true;

	// TODO unimplemented (never fails)

	// Non-planar case
	int match[3] = {-1, -1, -1};
	for (auto i : range(3))
		for (auto j : range(3))
			if (a[i] == b[j])
				match[i] = j;
	int count = 0; // how many vertices do two triangles have in common?
	for (auto i : range(3))
		if (match[i] >= 0)
			count += 1;
	if (count == 1) {

	}

	// Planar case
	// TODO there is some intersection:
	// OK case is if intersection is one vertex of both A and B (and no overlap)
	// OK case is if intersection is one edge of both A and B (and no overlap)
	return true;
}

Validity is_valid(const imesh3& mesh) {
	if (mesh.size() < 4)
		return Validity::TooFewFaces;

	// Minimal length of edge
	for (const auto& f : mesh)
		for (auto [a, b] : f.edges())
			if (a == b)
				return Validity::EdgeTooShort;

	// TODO any vertex can't be inside any other edge

	// Extract all edges
	std::unordered_map<isegment3, int> edge_to_face;
	for (auto i : range(mesh.size()))
		for (auto e : mesh[i].edges())
 			edge_to_face[e] = i;

	// Every edge must appear exactly twice (in opposite directions)
	for (auto [e, fi] : edge_to_face)
		if (edge_to_face.count(e.reversed()) == 0)
			return Validity::OpenEdge;

	// All triangles must form a single component
	std::vector<UnionFind> component(mesh.size());
	for (auto i : range(mesh.size()))
		for (auto e : mesh[i].edges())
			component[edge_to_face[e]].merge(component[i]);
	for (auto& c : component)
		if (c.find() != component[0].find())
			return Validity::SeparateComponents;

	// No self intersections
	for (auto i : range(mesh.size()))
		for (auto j : range(i))
			if (!are_valid_mesh_faces(mesh[i], mesh[j]))
				return Validity::SelfIntersection;

	return Validity::OK;
}

void make_valid(imesh3& mesh) {
	throw std::runtime_error("not implemented");
}
