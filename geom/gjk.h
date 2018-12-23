#include <geom/segment.h>

// 2D
bool AreConvexHullsIntersecting(cspan<double2> p, cspan<double2> q, double2* axis);
double DistanceBetweenConvexHulls(cspan<double2> p, cspan<double2> q, double2* axis);
double SignedDistanceBetweenConvexHulls(cspan<double2> p, cspan<double2> q, double2* axis);
optional<segment2> MinSegmentBetweenConvexHulls(cspan<double2> p, cspan<double2> q, double2* axis);

// 3D
bool AreConvexHullsIntersecting(cspan<double4> p, cspan<double4> q, double2* axis);
double DistanceBetweenConvexHulls(cspan<double4> p, cspan<double4> q, double2* axis);
double SignedDistanceBetweenConvexHulls(cspan<double4> p, cspan<double4> q, double2* axis);
optional<segment3> MinSegmentBetweenConvexHulls(cspan<double4> p, cspan<double4> q, double2* axis);
