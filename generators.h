#pragma once

#include "triangle.h"
#include "convex_hull.h"
#include "primitives.h"

imesh3 generate_box(int sx, int sy, int sz);

template<typename RND>
imesh3 generate_sphere(int vertices, int radius, RND& rnd) {
	// generate N random vertices on sphere
	std::vector<dvec3> V;
	V.resize(vertices);
	for (auto i : range(vertices))
		V[i] = random_direction(rnd);

	// increase the distance between the closest two vertices
	// (as random clumps vertices together)
	std::vector<dvec3> delta;
	delta.resize(vertices);
	for (auto e : range(40)) {
		for (dvec3& v : delta)
			v = {0, 0, 0};
		for (auto i : range(vertices))
			for (auto j : range(i)) {
				dvec3 d = V[i] - V[j];
				d *= 0.05 / glm::dot(d, d);
				delta[i] += d;
				delta[j] -= d;
			}
		for (auto i : range(vertices))
			V[i] = glm::normalize(V[i] + delta[i]);
	}

	// convert dvec3 to ivec3
	std::vector<ivec3> I;
	I.resize(vertices);
	for (auto i : range(vertices))
		I[i] = round(V[i] * (double)radius);

	return convex_hull(I);
}

imesh3 generate_cross(int inner, int outer);

// cone (rmin or rmax = 0) or cylinder
imesh3 generate_cylinder(int sides, int rmin, int rmax, int zmin, int zmax);

// platonic solids
imesh3 generate_tetrahedron(int radius);
imesh3 generate_cube(int radius);
imesh3 generate_octahedron(int radius);
imesh3 generate_dodecahedron(int radius);
imesh3 generate_icosahedron(int radius);

imesh3 generate_prism(const ipolygon2& poly, int zmin, int zmax);

imesh3 generate_pipe(const std::vector<ivec3>& path, int radius, int vertices);
