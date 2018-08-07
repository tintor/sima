#pragma once

#include "triangle.h"

// Volume of a valid polyhedron
double volume(const imesh3& mesh);
double signed_volume(const imesh3& mesh);
double triangle_volume(ivec3 a, ivec3 b, ivec3 c);

// Center of mass of a valid polyhedron
dvec3 center_of_mass(const imesh3& mesh);

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
//dmat3 moment_of_inertia(const imesh3& mesh);

bool is_aabb(const imesh3& mesh);
