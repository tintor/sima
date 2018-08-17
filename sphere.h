#pragma once
#include "format.h"
#include "span.h"

using sphere = double4;

inline sphere make_sphere(double3 center, double radius) {
	return {center.x, center.y, center.z, radius}; }
inline double3 center(sphere s) { return s.xyz; }
inline double radius(sphere s) { return s.w; }

sphere merge_spheres(sphere a, sphere b);
sphere bounding_sphere_minimal(span<const double3> points);
sphere bounding_sphere_fast(span<const double3> points);
