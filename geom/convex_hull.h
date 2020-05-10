#pragma once
#include <core/span.h>
#include <geom/mesh.h>

// Is this valid mesh a convex polyhedron?
bool is_convex(cspan<triangle3> mesh);

// Returns empty vector if no solution (points are coplanar)
void convex_hull(cspan<double3> points, mesh3& hull);

inline mesh3 convex_hull(cspan<double3> points) {
    mesh3 hull;
    convex_hull(points, /*out*/ hull);
    return hull;
}
