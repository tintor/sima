#include "sphere.h"
#include "vector.h"

sphere merge_spheres(sphere a, sphere b) {
	auto d = center(b - a);
	auto d2 = squared(d);

	if (squared(radius(a) - radius(b)) >= d2)
		return (radius(a) > radius(b)) ? a : b;

	auto dist = sqrt(d2);
	double c_radius = (radius(a) + radius(b) + dist) / 2;
	// division is safe here for small d, as center will converge to a.center
	double3 c_center = center(a) + ((c_radius - radius(a)) / dist) * d;
	return make_sphere(c_center, c_radius);
}
