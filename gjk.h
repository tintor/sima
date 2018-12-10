#include "segment.h"

bool AreConvexHullsIntersecting(span<const double2> ca, span<const double2> cb, double2* axis);
double DistanceBetweenConvexHulls(span<const double2> ca, span<const double2> cb, double2* axis);
double SignedDistanceBetweenConvexHulls(span<const double2> ca, span<const double2> cb, double2* axis);
optional<segment2> MinSegmentBetweenConvexHulls(span<const double2> ca, span<const double2> cb, double2* axis);
