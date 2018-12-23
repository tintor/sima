#pragma once
#include <geom/vector.h>
#include <geom/quaternion.h>
#include <core/exception.h>

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

struct Pose3 {
	double4 position; // w = 1
	quat orientation; // unit quaternion

	Pose3 inverse() {
		THROW(not_implemented);
	}

	// returned w will be 0
	double4 applyV(double4 v) const {
		return quat_rotate(orientation, v);
	}

	double4 apply(double4 p) const {
		return applyV(p) + position;
	}
};

// TODO combine poses
