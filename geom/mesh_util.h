#pragma once

#include <geom/triangle.h>

// find longest edge and two faces on it -> split the edge and the two triangles in half
// repeat until longest edge is short enough
void mesh_split_long_edges(imesh3& m, uint max_length);

// find and remove very flat tetrapacks (generalize for any triangle / polygon)
void mesh_reduce(imesh3& m, uint tolerance);

// round vertices more and more until two vertices merge, vertex merges with edge or vertex merges with face
void mesh_round(imesh3& a, uint amount);

// extrude a set of faces in some direction
void extrude_faces(imesh3& a, const vector<uint>& faces, ivec3 amount);
