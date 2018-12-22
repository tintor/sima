#pragma once
#include <core/format.h>
#include <core/span.h>
#include "vector.h"

class sphere {
public:
	sphere(double4 center, double radius) { s = center; s.w = radius; }
	point3 center() const { double4 e = s; e.w = 1; return e; }
	double radius() const { return s.w; }
	bool contains(point3 p) const { return squared(center() - p) <= squared(radius()); }
private:
	double4 s;
};

sphere minimal_sphere(sphere a, sphere b);
sphere minimal_sphere(sphere a, point3 b);
sphere minimal_sphere(point3 a, point3 b);
sphere minimal_sphere(point3 a, point3 b, point3 c);
sphere minimal_sphere(point3 a, point3 b, point3 c, point3 d);

sphere minimal_sphere(cspan<point3> points);

// not minimal, but close to it and faster to compute than minimal
sphere bounding_sphere(cspan<point3> points);
