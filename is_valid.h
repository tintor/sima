#pragma once
#include "triangle.h"
#include "format.h"

bool is_valid(const polygon2& poly);

enum class Validity {
	OK = 0,
	TooFewFaces = 1,
	InvalidFace = 2,
	OpenEdge = 3,
	SeparateComponents = 4,
	SelfIntersection = 5,
	NotSealed = 6,
	Inverted = 7,
	OverlappingFaces = 8,
};

Validity is_valid(const mesh3& mesh);

Validity is_valid(const mesh3& mesh);

void make_valid(mesh3& mesh);
