#pragma once
#include "triangle.h"
#include "span.h"

// Is this valid mesh a convex polyhedron?
bool is_convex(const mesh3& mesh);

// Returns empty vector if no solution (points are coplanar)
void convex_hull(span<const double3> points, mesh3& hull);

inline mesh3 convex_hull(span<const double3> points) {
	mesh3 hull;
	convex_hull(points, /*out*/hull);
	return hull;
}
