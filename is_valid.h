#pragma once
#include "mesh.h"
#include "triangle.h" // for polygon2
#include "format.h"

bool is_valid(const polygon2& poly);
bool is_valid(const xpolygon2& poly);

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
	NonPlanarFace = 9,
};

Validity is_valid(const mesh3& mesh);
Validity is_valid(const xmesh3& mesh);

void make_valid(mesh3& m);
void make_valid(xmesh3& m);
