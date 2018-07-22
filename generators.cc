#include "generators.h"

imesh3 generate_box(int sx, int sy, int sz) {
	std::vector<ivec3> vertices;
	vertices.reserve(8);
	for (int x : {-1, 1})
		for (int y : {-1, 1})
			for (int z : {-1, 1})
				vertices.emplace_back(x * sx, y * sy, z * sz);
	return convex_hull(vertices);
}

imesh3 generate_cylinder(int sides, int rmin, int rmax, int zmin, int zmax) {
	std::vector<ivec3> vertices;
	for (auto [z, r] : {std::pair{zmin, rmin}, std::pair{zmax, rmax}})
		if (r == 0)
			vertices.emplace_back(0, 0, z);
		else
			for (int i : range(sides)) {
				double a = (2 * M_PI / sides) * i;
				int x = std::round(std::cos(a) * r);
				int y = std::round(std::sin(a) * r);
				vertices.emplace_back(x, y, z);
			}
	return convex_hull(vertices);
}

void AddQuad(imesh3& mesh, ivec3 a, ivec3 b, ivec3 c, ivec3 d) {
	mesh.emplace_back(a, b, c);
	mesh.emplace_back(c, d, a);
}

imesh3 CreateCrossMesh(int inner, int outer) {
	imesh3 m;
	// TODO finish
	int a = inner, b = outer;

	// square face - Z axis
	AddQuad(m, ivec3(a, a, b), ivec3(-a, a, b), ivec3(-a, -a, b), ivec3(a, -a, b));
	AddQuad(m, ivec3(a, a, -b), ivec3(-a, -a, -b), ivec3(-a, a, -b), ivec3(a, -a, -b));
	// square face - Y axis
	AddQuad(m, ivec3(a, b, a), ivec3(-a, b, a), ivec3(-a, b, -a), ivec3(a, b, -a));
	AddQuad(m, ivec3(a, -b, a), ivec3(-a, -a, -a), ivec3(-a, -b, a), ivec3(a, -b, -a));
	// square face - X axis
	AddQuad(m, ivec3(b, a, a), ivec3(b, -a, a), ivec3(b, -a, -a), ivec3(b, a, -a));
	AddQuad(m, ivec3(-b, a, a), ivec3(-b, -a, -a), ivec3(-b, -a, a), ivec3(-b, a, -a));
	return m;
}
