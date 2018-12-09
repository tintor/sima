#pragma once
#include "polygon.h"
#include "mesh.h"

// TODO tessellation of 3d concave polygons with holes

// ear cutting algorithm
// assumes that P is valid and has not duplicate vertices or overlapping edges
void tesselate(polygon2 poly, /*out*/mesh2& tess);

inline mesh2 tesselate(polygon2 poly) {
	mesh2 tess;
	tesselate(poly, /*out*/tess);
	return tess;
}

// TODO circle cutting algorithm - to limit the maximum triangle size
