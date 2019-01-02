#pragma once
#include <geom/vector.h>
#include <geom/quaternion.h>
#include <core/exception.h>

struct Pose2 {
	double2 position;
	quat2 orientation; // unit vector

	Pose2() { }

	Pose2(double2 p, quat2 o) : position(p), orientation(o) { }

	// which inverse? there are two
	Pose2 inverse() {
		THROW(not_implemented);
	}

	double2 rotate(double2 v) const {
		return quat_rotate(orientation, v);
/*		return { orientation.x * v.x + orientation.y * v.y,
				-orientation.y * v.x + orientation.x * v.y};*/
	}

	double2 apply(double2 p) const {
		return rotate(p) + position;
	}
};

inline Pose2 interpolate(const Pose2& a, const Pose2& b, double t) {
	return Pose2(
		a.position + (b.position - a.position) * t,
		slerp(a.orientation, b.orientation, t));
}

inline Pose2 mat_to_pose(const mat33& a) {
	THROW(not_implemented);
}

inline mat33 pose_to_mat(const Pose2& a) {
	THROW(not_implemented);
}

inline mat33 mul(const mat33& p, const mat33& q) {
	mat33 e;
	e.a = mul(p, q.a);
	e.b = mul(p, q.b);
	e.c = mul(p, q.c);
	return e;
}

inline Pose2 mul(const Pose2& a, const Pose2& b) {
	return mat_to_pose(mul(pose_to_mat(a), pose_to_mat(b)));
}

struct Pose3 {
	double4 position; // w = 1
	quat orientation; // unit quaternion

	Pose3() { }

	Pose3(double4 p, quat o) : position(p), orientation(o) { }

	Pose3 inverse() {
		THROW(not_implemented);
	}

	// returned w will be 0
	double4 rotate(double4 v) const {
		double3 e = quat_rotate(orientation, v.xyz);
		return {e.x, e.y, e.z, 0};
	}

	double4 apply(double4 p) const {
		return rotate(p) + position;
	}
};

inline Pose3 interpolate(const Pose3& a, const Pose3& b, double t) {
	return Pose3(
		a.position + (b.position - a.position) * t,
		slerp(a.orientation, b.orientation, t));
}

inline Pose3 mat_to_pose(const mat44& a) {
	THROW(not_implemented);
}

inline mat44 pose_to_mat(const Pose3& a) {
	THROW(not_implemented);
}

inline mat44 mul(const mat44& p, const mat44& q) {
	mat44 e;
	e.a = mul(p, q.a);
	e.b = mul(p, q.b);
	e.c = mul(p, q.c);
	e.d = mul(p, q.d);
	return e;
}

inline Pose3 mul(const Pose3& a, const Pose3& b) {
	return mat_to_pose(mul(pose_to_mat(a), pose_to_mat(b)));
}
