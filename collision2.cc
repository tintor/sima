#include <core/dynamic_array.h>
#include <core/format.h>
#include <geom/classify.h>
#include <geom/gjk.h>
#include <geom/pose.h>
#include <geom/properties.h>
#include <geom/triangle.h>

struct ContactPoint2 {
    double2 position;
    double depth;
};

// Contact manifold describes the contact between two convex objects:
// Distinguishes between reference face and incident face:
// - reference face provides normal
// - incident face provides contact points
struct Contact2 {
    bool first_is_reference;  // iff the normal is pointing from first to second
    double2 normal;
    // TODO static_vector type
    ContactPoint2 points[2];
    int num_points;
    double time;  // relative
};

enum class ShapeType { Sphere, Capsule, Box, Polygon };

struct Shape2 {
    ShapeType type;
    double radius;
    double2 size;                     // for [axis aligned] box and capsule (rounded parts are on x axis)
    dynamic_array<double2> vertices;  // for polygon
                                      // TODO for polygon we also need list of planes with unit normals
};

template <typename T>
optional<pair<T, T>> solve_quadratic_equation(T a, T b, T c) {
    auto det = b * b - 4 * a * c;
    if (det < 0) return nullopt;
    auto s = sqrt(det);
    return {(-b - s) / (a * 2), (-b + s) / (a * 2)};
}

template <typename T>
optional<pair<T, T>> solve_quadratic_equation_no_div(T a, T b, T c) {
    auto det = b * b - 4 * a * c;
    if (det < 0) return nullopt;
    auto s = sqrt(det);
    return pair<T, T>{-b - s, -b + s};
}

// find first and last point on segment when it intersects circle centered in origin
template <typename Vec>
optional<pair<double, double>> segment_vs_circle(segment<Vec> s, double r, const double min_segment_squared = 1e-12) {
    Vec d = s.b - s.a;
    double dd = dot(d, d);
    if (dd < min_segment_squared) {
	if (dot(s.a, s.a) > r * r) return nullopt;
	return pair<double, double>{0, 1};
    }

    auto p = solve_quadratic_equation_no_div(dd, 2 * dot(s.a, d), dot(s.a, s.a) - r * r);
    double div = 2 * dd;
    if (!p || p->second < 0 || p->first > div) return nullopt;

    p->first = (p->first < 0) ? 0 : (p->first / div);
    p->second = (p->second > div) ? 1 : (p->second / div);
    return p;
}

// contact.time will be set to 1
bool shallow_collision(const Shape2& shape_a, const Shape2& shape_b, const pose2 pose_a, const pose2 pose_b,
                       Contact2& contact) {
    // quick bounding sphere test
    double r = shape_a.radius + shape_b.radius;
    if (squared(pose_a.position - pose_b.position) > r * r) return false;

    if (shape_a.type == ShapeType::Sphere && shape_b.type == ShapeType::Sphere) {
	// TODO contact normal and point
	double2 delta = pose_b.position - pose_a.position;
	double d = length(delta);
	double t = ((shape_a.radius - shape_b.radius) / d + 1) * 0.5;
	contact.points[0].position = pose_a.position * t + pose_b.position * (1 - t);
	contact.points[0].depth = 0;  // TODO finish
	contact.num_points = 1;
	contact.normal = delta / d;
	contact.time = 1;
	return true;
    }
    if (shape_a.type == ShapeType::Capsule && shape_b.type == ShapeType::Capsule) {
	// TODO contact normal and point
	contact.time = 1;
	return true;
    }

    THROW(not_implemented);
}

double2 support(double2 dir, const Shape2& shape) { THROW(not_implemented); }

double2 support(double2 dir, const Shape2& shape, pose2 pose) { THROW(not_implemented); }

// Assumes that objects rotate slightly between start and end poses
// Times are relative: start_a is at time 0, while end_a is at time 1
bool temporal_collision(const Shape2& shape_a, const Shape2& shape_b, pose2 start_a, pose2 start_b, pose2 end_a,
                        pose2 end_b, double start_time, double end_time, double max_depth) {
    double2 dir = {1, 0};  // TODO init from last frame
    return gjk_classify(
               [&shape_a, &shape_b, &start_a, &start_b, &end_a, &end_b, start_time, end_time](double2 d) {
	           // TODO buffer both objects by -max_depth
	           double2 start = support(d, shape_a, interpolate(start_a, end_a, start_time));
	           start -= support(-d, shape_b, interpolate(start_b, end_b, start_time));
	           double2 end = support(d, shape_a, interpolate(start_a, end_a, end_time));
	           end -= support(-d, shape_b, interpolate(start_b, end_b, end_time));
	           return (dot(d, start) > dot(d, end)) ? start : end;
               },
               [](double2 d) {
	           return double2{0, 0};
               },
               dir, 20) >= 0;
}

// Checks for deep collision at time [0, 1], and returns the earliest time
// Prevents tunneling of high speed objects.
bool deep_collision(const Shape2& shape_a, const Shape2& shape_b, pose2 start_a, pose2 start_b, pose2 end_a,
                    pose2 end_b, double max_depth, Contact2& contact) {
    // (1) check if two moving spheres are deeply colliding
    segment2 s(start_a.position - start_b.position, end_a.position - end_b.position);
    double r_deep = shape_a.radius + shape_b.radius - max_depth;
    // find min and max (relative) time when spheres are intersecting
    auto time = segment_vs_circle(s, r_deep);
    if (r_deep <= 0 || !time) return false;
    if (shape_a.type == ShapeType::Sphere && shape_b.type == ShapeType::Sphere) {
	// TODO contact normal and point
	contact.time = time->first;
	return true;
    }

    // (2) quick sweeping minkowski difference test (when no deep collision)
    double abs_error_a = shape_a.radius * (1 - cos(angle(start_a, end_a) / 2));
    double abs_error_b = shape_b.radius * (1 - cos(angle(start_b, end_b) / 2));
    double abs_error = abs_error_a + abs_error_b;

    // TODO what is acceptable error? and how to compute angle based on max error?
    const int n = 2;  // TODO

    double time_delta = (time->second - time->first) / n;
    for (auto i : range(n)) {
	double min_time = time->first + time_delta * i;
	double max_time = min_time + time_delta;

	if (temporal_collision(shape_a, shape_b, start_a, start_b, end_a, end_b, min_time, max_time, max_depth)) {
	    // (3) detailed sweeping minkowski difference test to find the deep collision time
	    while (max_time - min_time > 1e-3) {
		double mid_time = (min_time + max_time) / 2;
		if (temporal_collision(shape_a, shape_b, start_a, start_b, end_a, end_b, min_time, mid_time,
		                       max_depth)) {
		    max_time = mid_time;
		} else {
		    min_time = mid_time;
		}
	    }
	    // TODO contact normal and points
	    contact.time = (min_time + max_time) / 2;
	    return true;
	}
    }
    return false;
}
