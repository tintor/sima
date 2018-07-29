#include "generators.h"
#include "tesselate.h"

imesh3 generate_box(int sx, int sy, int sz) {
	std::vector<ivec3> vertices;
	vertices.reserve(8);
	for (int x : {-1, 1})
		for (int y : {-1, 1})
			for (int z : {-1, 1})
				vertices.push_back(ivec3{x * sx, y * sy, z * sz});
	return convex_hull(vertices);
}

imesh3 generate_cylinder(int sides, int rmin, int rmax, int zmin, int zmax) {
	std::vector<ivec3> vertices;
	for (auto [z, r] : {std::pair{zmin, rmin}, std::pair{zmax, rmax}})
		if (r == 0)
			vertices.push_back(ivec3{0, 0, z});
		else
			for (int i : range(sides)) {
				double a = (2 * M_PI / sides) * i;
				int x = std::round(std::cos(a) * r);
				int y = std::round(std::sin(a) * r);
				vertices.push_back(ivec3{x, y, z});
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
	AddQuad(m, ivec3{a, a, b}, ivec3{-a, a, b}, ivec3{-a, -a, b}, ivec3{a, -a, b});
	AddQuad(m, ivec3{a, a, -b}, ivec3{-a, -a, -b}, ivec3{-a, a, -b}, ivec3{a, -a, -b});
	// square face - Y axis
	AddQuad(m, ivec3{a, b, a}, ivec3{-a, b, a}, ivec3{-a, b, -a}, ivec3{a, b, -a});
	AddQuad(m, ivec3{a, -b, a}, ivec3{-a, -a, -a}, ivec3{-a, -b, a}, ivec3{a, -b, -a});
	// square face - X axis
	AddQuad(m, ivec3{b, a, a}, ivec3{b, -a, a}, ivec3{b, -a, -a}, ivec3{b, a, -a});
	AddQuad(m, ivec3{-b, a, a}, ivec3{-b, -a, -a}, ivec3{-b, -a, a}, ivec3{-b, a, -a});
	return m;
}

imesh3 generate_prism(const ipolygon2& poly, int zmin, int zmax) {
	imesh2 m2 = tesselate(poly);
	imesh3 m3;
	m3.reserve(m2.size() * 2 + poly.size() * 2);
	for (itriangle2 m : m2) {
		// TODO check orientation of m
		m3.emplace_back(ivec3{m.a.x, m.a.y, zmin}, ivec3{m.b.x, m.b.y, zmin}, ivec3{m.c.x, m.c.y, zmin});
		m3.emplace_back(ivec3{m.b.x, m.b.y, zmax}, ivec3{m.a.x, m.a.y, zmax}, ivec3{m.c.x, m.c.y, zmax});
	}
	ivec2 sa = poly.back();
	for (ivec2 sb : poly) {
		ivec3 a{sa.x, sa.y, zmin};
		ivec3 b{sa.x, sa.y, zmax};
		ivec3 c{sb.x, sb.y, zmax};
		ivec3 d{sb.x, sb.y, zmin};
		m3.emplace_back(a, b, c);
		m3.emplace_back(c, d, a);
		sa = sb;
	}
	return m3;
}

// TODO sphere with bumps: generate sphere, then move vertices randomly towards or away from center

// TODO concave imesh3 generator
// "grow" it from a starting shape, by either extruding a triangle (or two) or adding a pyramid over a face (or two),
// applied recursively (as long as there are no self intersections).
//
// "grow" a tree: first extrude face to get a branch, then add a pyramid to create a fork, then extrude again!
