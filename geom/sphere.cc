#include <geom/sphere.h>
#include <core/range.h>
#include <geom/aabb.h>
#include <core/exception.h>

sphere minimal_sphere(sphere a, sphere b) {
	double4 d = b.center() - a.center();
	double d2 = squared(d);

	if (squared(a.radius() - b.radius()) >= d2)
		return (a.radius() > b.radius()) ? a : b;

	auto dist = sqrt(d2);
	double c_radius = (a.radius() + a.radius() + dist) / 2;
	// division is safe here for small d, as center will converge to a.center
	point3 c_center = a.center() + ((c_radius - a.radius()) / dist) * d;
	return sphere(c_center, c_radius);
}

sphere minimal_sphere(sphere a, point3 b) {
	double4 d = a.center() - b;
	double d2 = squared(d);

	if (squared(a.radius()) >= d2)
		return a;

	auto dist = sqrt(d2);
	double c_radius = (a.radius() + dist) / 2;
	// division is safe here for small d, as center will converge to a.center
	point3 c_center = a.center() + ((c_radius - a.radius()) / dist) * d;
	return sphere(c_center, c_radius);
}

sphere correct_radius(sphere s, cspan<point3> points) {
	double r = s.radius();
	for (auto p : points)
		while (!s.contains(p)) {
			r = std::nexttoward(r, std::numeric_limits<double>::max());
			s = sphere(s.center(), r);
		}
	return s;
}

sphere sphere_from(point3 a, point3 b) {
	point3 center = avg(a, b);
	return sphere(center, length(a - center));
}

sphere minimal_sphere(point3 a, point3 b) {
	return correct_radius(sphere_from(a, b), {a, b});
}

sphere sphere_from(point3 a, point3 b, point3 c) {
    double4 A = a - c, B = b - c, X = cross(A, B);
    point3 center = c + cross(squared(A) * B - squared(B) * A, X) / (2 * squared(X));
	return sphere(center, length(a - center));
}

sphere minimal_sphere(point3 a, point3 b, point3 c) {
	sphere s = minimal_sphere(a, b);
	if (s.contains(c))
		return s;
	s = minimal_sphere(a, c);
	if (s.contains(b))
		return s;
	s = minimal_sphere(b, c);
	if (s.contains(a))
		return s;

	return correct_radius(sphere_from(a, b, c), {a, b, c});
}

// Finds E such that:
// x[i] * e.x + y[i] * e.y + z[i] * e.z = w[i]
double4 solve_linear_col(double4 x, double4 y, double4 z, double4 w) {
	double4 e = {det(w, y, z), det(x, w, z), det(x, y, w)};
	return e / det(x, y, z);
}

// Finds E such that:
// dot(a, e) = w[0]
// dot(b, e) = w[1]
// dot(c, e) = w[2]
double4 solve_linear_row(double4 a, double4 b, double4 c, double4 w) {
	double4 x = {a.x, b.x, c.x};
	double4 y = {a.y, b.y, c.y};
	double4 z = {a.z, b.z, c.z};
	return solve_linear_col(x, y, z, w);
}

sphere sphere_from(point3 a, point3 b, point3 c, point3 d) {
	// 2(x1 - x2) `dot` c = |x1|^2 - |x2|^2
	// 2(x1 - x3) `dot` c = |x1|^2 - |x3|^2
	// 2(x1 - x4) `dot` c = |x1|^2 - |x4|^2
	double sa = squared(a), sb = squared(b), sc = squared(c), sd = squared(d);
   	double4 w = {sa - sb, sa - sc, sa - sd};
	point3 center(solve_linear_row(2 * (a - b), 2 * (a - c), 2 * (a - d), w));
	center.w = 1;
	return sphere(center, length(a - center));
}

sphere minimal_sphere(point3 a, point3 b, point3 c, point3 d) {
	sphere s = minimal_sphere(a, b, c);
	if (s.contains(d))
		return s;
	s = minimal_sphere(a, b, d);
	if (s.contains(c))
		return s;
	s = minimal_sphere(a, c, d);
	if (s.contains(b))
		return s;
	s = minimal_sphere(b, c, d);
	if (s.contains(a))
		return s;

	return correct_radius(sphere_from(a, b, c, d), {a, b, c, d});
}

sphere minimal_sphere_brute_force_x(cspan<point3> points) {
	sphere minimal(double4{0, 0, 0, 1}, -1);
	for (auto a : range(points.size()))
		for (auto b : range(a + 1, points.size()))
			for (auto c : range(b + 1, points.size()))
				for (auto d : range(c + 1, points.size())) {
					sphere s = minimal_sphere(points[a], points[b], points[c], points[d]);
					if (s.radius() > minimal.radius())
						minimal = s;
				}
	return minimal;
}

sphere minimal_sphere(cspan<point3> interior, array<point3, 4>& surface, int surface_size) {
	if (surface_size == 4)
		return sphere_from(surface[0], surface[1], surface[2], surface[3]);

	if (surface_size == 3 && interior.size() == 0)
		return sphere_from(surface[0], surface[1], surface[2]);

	if (surface_size == 2 && interior.size() == 0)
		return sphere_from(surface[0], surface[1]);

	if (surface_size == 1 && interior.size() == 0)
		return sphere(surface[0], 0);

	if (surface_size == 1 && interior.size() == 1)
		return sphere_from(surface[0], interior[0]);

	if (surface_size == 0 && interior.size() == 1)
		return sphere(interior[0], 0);

	if (surface_size == 0 && interior.size() == 2)
		return sphere_from(interior[0], interior[1]);

	auto a = interior[0];
	interior = interior.pop_front();

	auto s = minimal_sphere(interior, surface, surface_size);
	if (!s.contains(a)) {
		surface[surface_size] = a;
		s = minimal_sphere(interior, surface, surface_size + 1);
	}
	return s;
}

sphere minimal_sphere(cspan<point3> points) {
	array<point3, 4> surface;
	return minimal_sphere(points, surface, 0);
}

sphere extermal_points_optimal_sphere(cspan<point3> points, cspan<float4> normals) {
	if (points.empty())
		THROW(invalid_argument, "points must be non-empty");
	if (points.size() <= normals.size() * 2)
		return minimal_sphere(points);

	aligned_vector<point3> extremal;
	// TODO maybe swap loops and compute two normals per point at once
	for (auto n : normals) {
		int imin = 0, imax = 0;
		// use float4 dot product intrinsic (as there is no AVX double4
		// dot product intrinsic and double precision is not needed)
		auto zmin = dot(n, points[0]), zmax = zmin;
		for (int i = 1; i < points.size(); i++) {
			auto z = dot(n, points[i]);
			if (z < zmin) {
				imin = i;
				zmin = z;
			}
			if (z > zmax) {
				imax = i;
				zmax = z;
			}
		}
		extremal.push_back(points[imin]);
		extremal.push_back(points[imax]);
	}
	remove_dups(extremal); // TODO will this make it slower or faster?

	sphere s = minimal_sphere(extremal);
	for (auto p : points)
		s = minimal_sphere(s, p);
	return s;
}

// TODO: generate even number of normals and evenly spaced (using sphere generator)
array<float4, 13> g_normals = {
	float4{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0},
	{1, 1, 0, 0}, {1, 0, 1, 0}, {0, 1, 1, 0}, {1, -1, 0, 0}, {1, 0, -1, 0}, {0, 1, -1, 0},
	{1, 1, 1, 0}, {1, 1, -1, 0}, {1, -1, 1, 0}, {-1, 1, 1, 0}
};

sphere bounding_sphere(cspan<point3> points) {
	return extermal_points_optimal_sphere(points, g_normals);
}
