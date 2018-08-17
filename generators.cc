#include "generators.h"
#include "tesselate.h"

mesh3 generate_box(double3 size) {
	vector<double3> vertices;
	vertices.reserve(8);
	for (double x : {-1, 1})
		for (double y : {-1, 1})
			for (double z : {-1, 1})
				vertices.push_back(double3{x, y, z} * size);
	return convex_hull(vertices);
}

mesh3 generate_cylinder(uint sides, double rmin, double rmax, double zmin, double zmax) {
	vector<double3> vertices;
	for (auto [z, r] : {pair{zmin, rmin}, pair{zmax, rmax}})
		if (r == 0)
			vertices.push_back(double3{0, 0, z});
		else
			for (int i : range(sides)) {
				double a = (2 * M_PI / sides) * i;
				double x = std::cos(a) * r;
				double y = std::sin(a) * r;
				vertices.push_back(double3{x, y, z});
			}
	return convex_hull(vertices);
}

void AddQuad(mesh3& mesh, double3 a, double3 b, double3 c, double3 d) {
	mesh.emplace_back(a, b, c);
	mesh.emplace_back(c, d, a);
}

mesh3 CreateCrossMesh(double inner, double outer) {
	mesh3 m;
	// TODO finish
	double a = inner, b = outer;

	// square face - Z axis
	AddQuad(m, double3{a, a, b}, double3{-a, a, b}, double3{-a, -a, b}, double3{a, -a, b});
	AddQuad(m, double3{a, a, -b}, double3{-a, -a, -b}, double3{-a, a, -b}, double3{a, -a, -b});
	// square face - Y axis
	AddQuad(m, double3{a, b, a}, double3{-a, b, a}, double3{-a, b, -a}, double3{a, b, -a});
	AddQuad(m, double3{a, -b, a}, double3{-a, -a, -a}, double3{-a, -b, a}, double3{a, -b, -a});
	// square face - X axis
	AddQuad(m, double3{b, a, a}, double3{b, -a, a}, double3{b, -a, -a}, double3{b, a, -a});
	AddQuad(m, double3{-b, a, a}, double3{-b, -a, -a}, double3{-b, -a, a}, double3{-b, a, -a});
	return m;
}

mesh3 generate_prism(const polygon2& poly, double zmin, double zmax) {
	mesh2 m2 = tesselate(poly);
	mesh3 m3;
	m3.reserve(m2.size() * 2 + poly.size() * 2);
	for (triangle2 m : m2) {
		// TODO check orientation of m
		m3.emplace_back(double3{m.a.x, m.a.y, zmin}, double3{m.b.x, m.b.y, zmin}, double3{m.c.x, m.c.y, zmin});
		m3.emplace_back(double3{m.b.x, m.b.y, zmax}, double3{m.a.x, m.a.y, zmax}, double3{m.c.x, m.c.y, zmax});
	}
	double2 sa = poly.back();
	for (double2 sb : poly) {
		double3 a{sa.x, sa.y, zmin};
		double3 b{sa.x, sa.y, zmax};
		double3 c{sb.x, sb.y, zmax};
		double3 d{sb.x, sb.y, zmin};
		m3.emplace_back(a, b, c);
		m3.emplace_back(c, d, a);
		sa = sb;
	}
	return m3;
}

// TODO sphere with bumps: generate sphere, then move vertices randomly towards or away from center

// TODO concave mesh3 generator
// "grow" it from a starting shape, by either extruding a triangle (or two) or adding a pyramid over a face (or two),
// applied recursively (as long as there are no self intersections).
//
// "grow" a tree: first extrude face to get a branch, then add a pyramid to create a fork, then extrude again!
