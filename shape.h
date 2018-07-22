#pragma once

#include "rendering.h"
#include "primitives.h"
#include "triangle.h"
#include "solid_bsp_tree.h"
#include "transform.h"
#include "aabb.h"

// Shape is 3d solid, immutable, purely geometric and with origin in center of mass
struct Shape {
	Shape(const imesh3& mesh, /*out*/transform3& pose);

	const imesh3 faces;
	const std::vector<isegment3> convex_edges;
	const std::vector<ivec3> convex_vertices;
	
	const SolidBSPTree solid_bsp_tree;

	const dmat3 inertia_tensor; // assuming density of 1kg/m^3
	const double volume;
	const double sphere_radius;
	const aabb<ivec3> box;
	const bool is_convex;
	const bool is_box;
};

real signed_distance(const dvec3& v, const Shape& shape);

// either O(n*logn) SolidLeafBSPTree or O(n) algorithm using nearest face
inline bool intersects(const dvec3& v, const Shape& shape) {
	// TODO bounding box / sphere check
	// TODO optimize for box
	return signed_distance(v, shape) < 0;
}

real squared_distance_segment_origin(const segment3& p);

// TODO find all faces within epsilon distance from the edge
//      if any triangle is penetrated internally return true
//      if
//		rare_case:
//      collect at most one intersection point of edge with every nearby triangle
//		sort these points and check every mid point against the interior of shape
bool intersects_edge_interior(const segment3& edge, const Shape& shape);

bool intersects(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b);

struct Contact {
	real squared_dist;
	dvec3 position; // mid point
	dvec3 normal; // unit normal from body B to body A
};

bool is_vertex_triangle_contact(const dvec3& p, const triangle3& m, Contact& /*out*/contact);
bool is_edge_edge_contact(const segment3& p, const segment3& q, Contact& /*out*/contact);
// Assuming shapes aren't intersecting, look for pairs of features that are within ContactEpsilon
std::vector<Contact> find_all_contacts(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b);

// returns 0 if interecting
real distance(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b);
