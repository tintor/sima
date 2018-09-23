#include "classify.h"
#include "exception.h"
#include "edges.h"
#include "project.h"
#include "properties.h"
#include "sphere.h"
#include "primitives.h"

double Distance(sphere s, double4 v) {
	return max(0.0, length(s.center() - v) - s.radius());
}
double Distance(sphere a, sphere b) {
	return max(0.0, length(a.center() - b.center()) - a.radius() - b.radius());
}

int ClassifyFacePrismPoint(const face& f, double4 v, plane p) {
	// TODO project to 2d, and run Classify
	THROW(not_implemented);
}

vector<segment3> UniqueEdges(const xmesh3& ma) {
	THROW(not_implemented);
}

inline segment3 ReverseIf(segment3 s, bool reverse) {
	return reverse ? s.reversed() : s;
}

static double DistanceFaceVertexPairs(const xmesh3& ma, const xmesh3& mb, sphere sb, double min_distance, bool reverse, segment3* nearest) {
	for (const face& f : ma) {
		optional<plane> p;
		auto s = bounding_sphere(f.vertices());

		if (Distance(s, sb) < min_distance)
			for (const double4& v : mb.vertices()) // TODO faster if iterating over unique vertices
				if (squared(s.center() - v) < squared(min_distance + s.radius())) {
					if (!p.has_value())
						p = best_fit_plane(f.vertices());
					if (abs(p->distance(v)) < min_distance && ClassifyFacePrismPoint(f, v, *p) <= 0) {
						double signed_distance = p->distance(v);
						double d = abs(signed_distance);
						if (d < min_distance) {
							if (nearest)
								*nearest = ReverseIf(segment3(v - signed_distance * p->normal(), v), reverse);
							min_distance = d;
						}
					}
				}
	}
	return min_distance;
}

// distance between the closest pair of points
double Distance(span<const double4> ma, span<const double4> mb, sphere sa, sphere sb, segment3* nearest) {
	// TODO can make it faster if there is a separating plane between two point sets
	//      use two bounding spheres to estimate potentiall separating plane
	//      sort both sets based on distance from plane
	//      start comparing nearest pairs and maintain lower_bound and upper_bound until remaining vertices are far from plane

	double ds = std::numeric_limits<double>::max();
	for (double4 a : ma)
		if (squared(a - sb.center()) < squared(ds + sb.radius()))
			for (double4 b : mb) {
				double e = squared(a - b);
				if (e < ds) {
					ds = e;
					if (nearest)
						*nearest = segment3(a, b);
				}
			}
	return sqrt(ds);
}

double Distance(const xmesh3& ma, const xmesh3& mb, segment3* nearest) {
	assert(Classify(ma, mb) >= 0);

	auto sa = bounding_sphere(ma.vertices());
	auto sb = bounding_sphere(mb.vertices());
	double min_distance = Distance(ma.vertices(), mb.vertices(), sa, sb, nearest);

	// distance between all edge pairs
	// TODO faster: swap ma / mb if sa is smaller
	for (const auto& ea : UniqueEdges(ma)) {
		auto s = minimal_sphere(ea.a, ea.b);
		if (squared(s.center() - sb.center()) < squared(min_distance + s.radius() + sb.radius()))
			for (const auto& eb : UniqueEdges(mb)) {
				auto w = ::nearest(ea, eb);
				double d = length(w.first.a - w.first.b);
				if (d < min_distance) {
					min_distance = d;
					if (nearest)
						*nearest = w.first;
				}
			}
	}

	min_distance = DistanceFaceVertexPairs(ma, mb, sb, min_distance, false, nearest);
	min_distance = DistanceFaceVertexPairs(mb, ma, sa, min_distance, true, nearest);
	return min_distance;
}
