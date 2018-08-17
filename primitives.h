#pragma once
#include "vector.h"
#include "segment.h"
#include "triangle.h"
#include "plane.h"

inline constexpr long double operator "" _deg(long double a) { return a * (M_PI / 180); }
inline constexpr long double operator "" _deg(unsigned long long a) { return a * (M_PI / 180); }

inline string deg(long double r) { return format("%g deg", r * (180 / M_PI)); }

inline constexpr double clamp(double t, double min = 0, double max = 1) {
	if (t < min) return min;
	if (t > max) return max;
	return t;
}

// How close a point needs to be to a plane to be considered on the plane?
constexpr double PlanarEpsilon = 1e-6;

// How close two features need to be to be considered touching?
// TODO can it be same as PlanarEpsilon?
constexpr double ContactEpsilon = 10 * PlanarEpsilon;

double squared_distance(segment3 p, segment3 q);
double squared_distance(line3 m, double3 p);

double distance(double3 a, double3 b);
double distance(double3 a, segment3 b);
double distance(segment3 a, double3 b);
double distance(segment3 a, segment3 b);

constexpr bool inside_triangle_prism(double3 p, triangle3 m, double3 normal);

double distance(double3 p, triangle3 m);
double distance(double3 v, triangle3 m, plane p_hint);
double distance(triangle3 m, double3 v);

bool intersects(line3 e, triangle3 m);
bool intersects2(line3 e, triangle3 m);
bool intersection(line3 e, triangle3 m, /*out*/double3& result);
bool intersects_in_point(segment3 e, triangle3 m);
bool intersects_in_point(ray3 e, triangle3 m);

double disjoint_distance(segment3 e, triangle3 m);
double distance(segment3 e, triangle3 m);

double distance(triangle3 m, segment3 e);
double disjoint_distance(triangle3 p, triangle3 q);
double distance(triangle3 p, triangle3 q);

// Angle between oriented triangles ABC and BAD.
// If edge is convex then angle will be <PI
// If edge is planar then angle will be =PI
// If edge is concave than angle will be >PI
double edge_angle(double3 a, double3 b, double3 c, double3 d);
