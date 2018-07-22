#pragma once

#include "triangle.h"

int128 signed_volume_mul6(const imesh3& mesh);

// Volume of a valid polyhedron
double volume(const imesh3& mesh);

// Center of mass of a valid polyhedron
ivec3 center_of_mass(const imesh3& mesh);

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
dmat3 moment_of_inertia(const imesh3& mesh);

bool is_aabb(const imesh3& mesh);
