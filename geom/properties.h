#pragma once
#include <geom/polygon.h>
#include <geom/plane.h>
#include <geom/mesh.h>

double SignedTriangleVolume6(double3 a, double3 b, double3 c);

double SignedVolume(const mesh3& mesh);
double SignedVolume(const xmesh3& mesh);

inline double Volume(const mesh3& mesh) { return std::abs(SignedVolume(mesh)); }
inline double Volume(const xmesh3& mesh) { return std::abs(SignedVolume(mesh)); }

double3 CenterOfMass(const mesh3& mesh);
double3 CenterOfMass(const xmesh3& mesh);

double2 CenterOfMass(const polygon2& poly);
double2 centroid(cspan<double2> poly);

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
//dmat3 moment_of_inertia(const imesh3& mesh);

bool is_aabb(const mesh3& mesh);

double3 eigen_vector(cspan<double3> points);

inline plane best_fit_plane(cspan<double3> points) {
	if (points.size() == 3) {
		return plane(points[0], points[1], points[2]);
	}
	double3 normal = normalize(eigen_vector(points));
	double s = 0;
	for (double3 p : points)
		s += dot(normal, p);
	return plane(normal, s / points.size());
}
