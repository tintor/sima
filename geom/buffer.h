#pragma once
#include <core/std.h>
#include <geom/vector.h>

// convex.size() > 0
// if convex.size() >= 3 it must be oriented counter-clockwise
// radius > 0
// circle_vertices >= 4
vector<double2> ComputeBuffer(cspan<double2> convex, double radius, int circle_vertices);
