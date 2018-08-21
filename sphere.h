#pragma once
#include "format.h"
#include "span.h"
#include "vector.h"

using sphere = double4;

inline sphere make_sphere(point3 center, double radius) { return {center.x(), center.y(), center.z(), radius}; }
inline point3 center(sphere s) { return s.xyz; }
inline double radius(sphere s) { return s.w; }
inline bool contains(sphere s, point3 p) { return squared(center(s) - p) <= squared(radius(s)); }

sphere minimal_sphere(sphere a, sphere b);
sphere minimal_sphere(sphere a, point3 b);
sphere minimal_sphere(point3 a, point3 b);
sphere minimal_sphere(point3 a, point3 b, point3 c);
sphere minimal_sphere(point3 a, point3 b, point3 c, point3 d);

sphere minimal_sphere(span<const point3> points);

// not minimal, but close to it and faster to compute than minimal
sphere bounding_sphere(span<const point3> points);
