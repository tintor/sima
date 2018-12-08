#pragma once
#include "segment.h"
#include "mesh.h"
#include "polygon.h"

// TODO project should preserve distances!!!
double2 Project(double4 v, int axis);
xpolygon2 Project(const face& f, int axis);
ray2 Project(const ray3& s, int axis);
segment2 Project(const segment3& s, int axis);
int ProjectionAxis(const face& f);
