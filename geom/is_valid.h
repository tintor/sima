#pragma once
#include <core/format.h>
#include <geom/mesh.h>
#include <geom/polygon.h>

bool IsValidPolygon(cspan<double2> poly);
bool IsValid(const polygon2& poly);
bool IsValid(const xpolygon2& poly);

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

Validity IsValid(const mesh3& mesh);
Validity IsValid(const xmesh3& mesh);

void MakeValid(mesh3& m);
void MakeValid(xmesh3& m);

bool AreAllEdgesConnected(const mesh3& mesh);
