#pragma once
#include "polygon.h"
#include "plane.h"
#include "mesh.h"

double SignedTriangleVolume6(double4 a, double4 b, double4 c);

double SignedVolume(const mesh3& mesh);
double SignedVolume(const xmesh3& mesh);

inline double Volume(const mesh3& mesh) { return std::abs(SignedVolume(mesh)); }
inline double Volume(const xmesh3& mesh) { return std::abs(SignedVolume(mesh)); }

double4 CenterOfMass(const mesh3& mesh);
double4 CenterOfMass(const xmesh3& mesh);

double2 CenterOfMass(const polygon2& poly);

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
//dmat3 moment_of_inertia(const imesh3& mesh);

bool is_aabb(const mesh3& mesh);

double4 eigen_vector(span<const point3> points);

inline plane best_fit_plane(span<const point3> points) {
	if (points.size() == 3) {
		return plane(points[0], points[1], points[2]);
	}
	double4 normal = normalize(eigen_vector(points));
	double s = 0;
	for (point3 p : points)
		s += dot(normal, p);
	return plane(normal, s / points.size());
}
