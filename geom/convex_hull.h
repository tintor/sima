#pragma once
#include <geom/mesh.h>
#include <core/span.h>

// Is this valid mesh a convex polyhedron?
bool is_convex(cspan<triangle3> mesh);

// Returns empty vector if no solution (points are coplanar)
void convex_hull(cspan<double4> points, mesh3& hull);

inline mesh3 convex_hull(cspan<double4> points) {
	mesh3 hull;
	convex_hull(points, /*out*/hull);
	return hull;
}
