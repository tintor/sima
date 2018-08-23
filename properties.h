#pragma once

#include "triangle.h"

// Volume of a valid polyhedron
double volume(const mesh3& mesh);
double signed_volume(const mesh3& mesh);
double triangle_volume(double3 a, double3 b, double3 c);

// Center of mass of a valid polyhedron
double3 center_of_mass(const mesh3& mesh);

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
//dmat3 moment_of_inertia(const imesh3& mesh);

bool is_aabb(const mesh3& mesh);

double3 eigen_vector(span<const double3> points);
