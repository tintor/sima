#pragma once
#include "triangle.h"

// ear cutting algorithm
// assumes that P is valid and has not duplicate vertices or overlapping edges
void tesselate(ipolygon2 poly, /*out*/imesh2& tess);

inline imesh2 tesselate(ipolygon2 poly) {
	imesh2 tess;
	tesselate(poly, /*out*/tess);
	return tess;
}

// TODO circle cutting algorithm - to limit the maximum triangle size
