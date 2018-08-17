#pragma once
#include "format.h"

struct plucker {
	double3 u, v;

	static plucker from_line(double3 a, double3 b) {
		return plucker{b - a, cross(b, a)};
	}

	// <0 Clockwise (if you look in direction of one line, other will go CW around it)
 	// =0 Intersect or Parallel
  	// >0 Counterclockwise
	double crossing(plucker p) const {
		return dot(u, p.v) + dot(v, p.u);
	}
};
