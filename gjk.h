#include "segment.h"

// 2D
bool AreConvexHullsIntersecting(span<const double2> p, span<const double2> q, double2* axis);
double DistanceBetweenConvexHulls(span<const double2> p, span<const double2> q, double2* axis);
double SignedDistanceBetweenConvexHulls(span<const double2> p, span<const double2> q, double2* axis);
optional<segment2> MinSegmentBetweenConvexHulls(span<const double2> p, span<const double2> q, double2* axis);

// 3D
bool AreConvexHullsIntersecting(span<const double4> p, span<const double4> q, double2* axis);
double DistanceBetweenConvexHulls(span<const double4> p, span<const double4> q, double2* axis);
double SignedDistanceBetweenConvexHulls(span<const double4> p, span<const double4> q, double2* axis);
optional<segment3> MinSegmentBetweenConvexHulls(span<const double4> p, span<const double4> q, double2* axis);
