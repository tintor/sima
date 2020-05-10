#pragma once
#include <core/hash.h>
#include <geom/matrix.h>
#include <geom/quaternion.h>

inline double2 rotate(double angle, double2 v) {
    double c = cos(angle);
    double s = sin(angle);
    return {c * v.x + s * v.y, -s * v.x + c * v.y};
}

struct pose2 {
    double2 position;
    double orientation;

    pose2() {}
    pose2(double2 p, double o) : position(p), orientation(o) {}

    double2 rotate(double2 v) const { return ::rotate(orientation, v); }
    double2 apply(double2 p) const { return rotate(p) + position; }

    bool operator==(pose2 o) const { return equal(position, o.position) && orientation == o.orientation; }
    bool operator!=(pose2 o) const { return !operator==(o); }
};

inline hash operator<<(hash h, pose2 p) { return h << p.position << p.orientation; }

inline double slerp(double angle_a, double angle_b, double t) {
    double d = angle_b - angle_a;
    while (d >= PI) d -= 2 * PI;
    while (d < -PI) d += 2 * PI;
    return angle_a + d * t;
}

inline pose2 interpolate(pose2 a, pose2 b, double t) {
    return pose2(a.position + (b.position - a.position) * t, slerp(a.orientation, b.orientation, t));
}

inline pose2 mat_to_pose(double33 m) {
    double22 e = {m.a.xy, m.b.xy};
    return pose2(m.c.xy, atan2(m.a.y, m.a.x));
}

inline double33 pose_to_mat(pose2 p) {
    double c = cos(p.orientation);
    double s = sin(p.orientation);
    return {{c, s, 0}, {-s, c, 0}, extend(p.position, 1)};
}

inline pose2 mul(pose2 a, pose2 b) { return pose2(b.apply(a.position), a.orientation + b.orientation); }

inline pose2 inv(pose2 a) { return pose2(-rotate(-a.orientation, a.position), -a.orientation); }

inline double angle(pose2 a, pose2 b) {
    double d = b.orientation - a.orientation;
    while (d >= PI) d -= 2 * PI;
    while (d < -PI) d += 2 * PI;
    return abs(d);
}

struct pose3 {
    double3 position;
    quat orientation;  // unit quaternion

    pose3() {}
    pose3(double3 p, quat o) : position(p), orientation(o) {}

    double3 rotate(double3 v) const { return quat_rotate(orientation, v.xyz); }
    double3 apply(double3 p) const { return rotate(p) + position; }

    bool operator==(pose3 o) const { return equal(position, o.position) && equal(orientation, o.orientation); }
    bool operator!=(pose3 o) const { return !operator==(o); }
};

inline hash operator<<(hash h, pose3 p) { return h << p.position << p.orientation; }

inline pose3 interpolate(pose3 a, pose3 b, double t) {
    return pose3(a.position + (b.position - a.position) * t, slerp(a.orientation, b.orientation, t));
}

inline pose3 mat_to_pose(const double44& m) {
    double33 e = {m.a.xyz, m.b.xyz, m.c.xyz};
    return pose3(m.d.xyz, quat_from_matrix(e));
}

inline double44 pose_to_mat(pose3 p) {
    double33 e = quat_to_matrix(p.orientation);
    return {extend(e.a, 0), extend(e.b, 0), extend(e.c, 0), extend(p.position, 1)};
}

inline pose3 mul(pose3 a, pose3 b) { return mat_to_pose(pose_to_mat(a) * pose_to_mat(b)); }

inline pose3 inverse(pose3 a) { return mat_to_pose(inv(pose_to_mat(a))); }

inline double angle(pose3 a, pose3 b) { return 2 * acos(dot(a.orientation, b.orientation)); }
