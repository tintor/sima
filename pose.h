#pragma once
#include "vector.h"
#include "exception.h"

struct Pose2 {
	double2 position;
	double2 orientation; // unit vector

	Pose2 inverse() {
		THROW(not_implemented);
	}

	double2 applyV(double2 v) const {
		return { orientation.x * v.x + orientation.y * v.y,
				-orientation.y * v.x + orientation.x * v.y};
	}

	double2 apply(double2 p) const {
		return applyV(p) + position;
	}
};

// TODO combine poses

inline double4 quat_mul_vec(double4 q, double4 v) {
	double x2 = q.x * 2;
	double y2 = q.y * 2;
	double z2 = q.z * 2;
	double xx2 = q.x * x2;
	double yy2 = q.y * y2;
	double zz2 = q.z * z2;
	double num7 = q.x * y2;
	double num8 = q.x * z2;
	double num9 = q.y * z2;
	double num10 = q.w * x2;
	double num11 = q.w * y2;
	double zw2 = q.w * z2;

	double4 mx = {1 - (yy2 + zz2), num7 + zw2, num8 - num11, 0};
	double4 my = {num7 - zw2, 1 - (xx2 + zz2), num9 + num10, 0};
	double4 mz = {num8 + num11, num9 - num10, 1 - (xx2 + yy2), 0};

	return mx * v.x + my * v.y + mz * v.z;
}

struct Pose3 {
	double4 position; // w = 1
	double4 orientation; // unit quaternion

	Pose3 inverse() {
		THROW(not_implemented);
	}

	// returned w will be 0
	double4 applyV(double4 v) const {
		return quat_mul_vec(orientation, v);
	}

	double4 apply(double4 p) const {
		return applyV(p) + position;
	}
};

// TODO combine poses
