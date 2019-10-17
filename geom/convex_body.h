#pragma once
#include <core/std.h>
#include <geom/mesh.h>

struct xyz4 {
	double4 x, y, z;
};

// mesh3 with extra precomputed fields for faster ClassifyConvexConvex!
struct cmesh3 {
	mesh3 mesh;
	vector<double3> vertices;
	vector<xyz4> vertices4;
	vector<double3> faceAxis;
	vector<double3> edgeAxis;
};

cmesh3 GenerateConvexMesh(const vector<double3>& vertices);

// in case of penetration return axis with min
// returns:
// 1 if disjoint
// 0 if in contact
// - returns axis of contact (oriented from ca to cb)
// - contacts represent convex 2D shape in 3D space where meshes are touching
// - vertex/face contact -> 1 contact point
// - edge/edge contact -> 1 contact point
// - edge/face contact -> 2 contact points
// - face/face contact -> 3+ contact points
// -1 if overlaping
// - returns axis of minimal overlap
// - returns amount of minimal overlap
int ClassifyConvexConvex(
	const cmesh3& ca,
	const cmesh3& cb,
	bool hasNormal,
	double3& normal,
	vector<double3>& contacts,
	double& overlap,
	vector<double3>& workA,
	vector<double3>& workB);
