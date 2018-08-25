#pragma once
#include "triangle.h"
#include "plane.h"

// Volume of a valid polyhedron
double volume(const mesh3& mesh);
double signed_volume(const mesh3& mesh);
double triangle_volume(double4 a, double4 b, double4 c);

// Center of mass of a valid polyhedron
double4 center_of_mass(const mesh3& mesh);

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
//dmat3 moment_of_inertia(const imesh3& mesh);

bool is_aabb(const mesh3& mesh);

double4 eigen_vector(span<const point3> points);

inline plane best_fit_plane(span<const point3> points) {
	double4 normal = normalize(eigen_vector(points));
	double s = 0;
	for (point3 p : points)
		s += dot(normal, p);
	return plane(normal, s / points.size());
}
