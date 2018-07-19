#pragma once

#include "triangle.h"

enum class Validity {
	OK = 0,
	EdgeTooShort = 1,
	TooFewFaces = 2,
	OpenEdge = 3,
	SeparateComponents = 4,
	SelfIntersection = 5,
};

Validity is_valid(const imesh3& mesh);

void make_valid(imesh3& mesh);
