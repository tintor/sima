#pragma once
#include "triangle.h"
#include "array_ptr.h"

// Is this valid mesh a convex polyhedron?
bool is_convex(const imesh3& mesh);

// Returns empty vector if no solution (points are coplanar)
void convex_hull(array_cptr<ivec3> points, imesh3& hull);

inline imesh3 convex_hull(array_cptr<ivec3> points) {
	imesh3 hull;
	convex_hull(points, /*out*/hull);
	return hull;
}
