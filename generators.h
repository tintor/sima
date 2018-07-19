#pragma once

#include "triangle.h"
#include "convex_hull.h"
#include "primitives.h"

imesh3 generate_box(int sx, int sy, int sz);

template<typename RND>
imesh3 generate_sphere(int vertices, int radius, RND& rnd) {
	std::vector<ivec3> V(vertices);
	for (auto i : range(vertices))
		V[i] = random_direction(rnd) * radius;
	return convex_hull(V);
}

imesh3 generate_cross(int inner, int outer);

// cone (rmin or rmax = 0) or cylinder
imesh3 generate_cylinder(int vertices, int rmin, int rmax, int zmin, int zmax);

// platonic solids
imesh3 generate_tetrahedron(int radius);
imesh3 generate_cube(int radius);
imesh3 generate_octahedron(int radius);
imesh3 generate_dodecahedron(int radius);
imesh3 generate_icosahedron(int radius);

imesh3 generate_prism(const ipolygon2& poly, int zmin, int zmax);

imesh3 generate_pipe(const std::vector<ivec3>& path, int radius, int vertices);
