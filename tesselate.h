#pragma once
#include "triangle.h"

// ear cutting algorithm
// assumes that P is valid and has not duplicate vertices or overlapping edges
imesh2 tesselate(ipolygon2 p);

// TODO circle cutting algorithm - to limit the maximum triangle size
