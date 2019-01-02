#pragma once
#include <geom/vector.h>
#include <core/exception.h>

using quat2 = double2;

inline quat2 quat_from_angle(double angle) {
	return {sin(angle / 2), cos(angle / 2)};
}

inline quat2 quat_mul(quat2 a, quat2 b) {
	return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
}

// returns q such that a * q == q * a == Identity
// (assumes that q is unit quaternion, so no need to normalize)
inline quat2 quat_inv(quat2 a) {
	return {-a.x, a.y};
}

// TODO is there a loss of precision due to acos?
inline double quat_angle(quat2 q) {
	return acos(q.y) * 2;
}

// single both are assumed to be unit quats, this function returns
// cos half angle to rotate from a to b
inline double quat_dot(quat2 a, quat2 b) {
	return a.y * b.y + a.x * b.x;
}

inline double2 quat_rotate(quat2 q, double2 a) {
	THROW(not_implemented);
}

// Linear Interpolation
inline quat2 lerp(quat2 p, quat2 q, double t) {
	return p * (1 - t) + q * t;
}

// Spherical Linear Interpolation
inline quat2 slerp(quat2 p, quat2 q, double t) {
	double d = dot(p, q);
	if (d >= 0) {
		double a = acos(d);
		double k = 1 / sin(a);
		return p * (sin(a - t * a) * k) + q * (sin(t * a) * k);
	}
	double a = acos(-d);
	double k = 1 / sin(a);
	return p * (sin(a - t * a) * k) - q * (sin(t * a) * k);
}

// Partialy based on Matrix and Quaternion FAQ, http://mccammon.ucsd.edu/~adcock/matrixfaq.html

using quat = double4;

constexpr quat IDENTITY_QUAT = {0, 0, 0, 1};

// Creates quaternion that represents rotation around axis by angle (in radians).
//  Right-hand rule is used for rotation direction.
inline quat quat_from_axis_angle(double4 axis, double angle) {
	auto q = axis * (sin(angle / 2) / length(axis));
	q.w = cos(angle / 2);
	return q;
}

// Creates quaternion that rotates vector a to vector b
inline quat quat_from_to(double4 a, double4 b) {
	auto axis = cross(a, b);
	if (squared(axis) < 1e-12)
		axis = any_normal(a);
	return quat_from_axis_angle(axis, acos(dot(a, b) / sqrt(squared(a) * squared(b))));
}

inline quat quat_from_x_axis(double angle) {
	return {sin(angle / 2), 0, 0, cos(angle / 2)};
}

inline quat quat_from_y_axis(double angle) {
	return {0, sin(angle / 2), 0, cos(angle / 2)};
}

inline quat quat_from_z_axis(double angle) {
	return {0, 0, sin(angle / 2), cos(angle / 2)};
}

inline quat quat_mul(quat a, quat b) {
	double x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
	double y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
	double z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
	double w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
	return {x, y, z, w};
}

// returns q such that a * q == q * a == Identity
// (assumes that q is unit quaternion, so no need to normalize)
inline quat quat_inv(quat a) {
	return {-a.x, -a.y, -a.z, a.w};
}

// inv(a) * b
inline double4 quat_ldiv(quat a, quat b) {
	double x = a.w * b.x - a.x * b.w - a.y * b.z + a.z * b.y;
	double y = a.w * b.y + a.x * b.z - a.y * b.w - a.z * b.x;
	double z = a.w * b.z - a.x * b.y + a.y * b.x - a.z * b.w;
	double w = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
	return {x, y, z, w};
}

// a * inv(b)
inline quat quat_rdiv(quat a, quat b) {
	double x = -a.w * b.x + a.x * b.w - a.y * b.z + a.z * b.y;
	double y = -a.w * b.y + a.x * b.z + a.y * b.w - a.z * b.x;
	double z = -a.w * b.z - a.x * b.y + a.y * b.x + a.z * b.w;
	double w = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
	return {x, y, z, w};
}

// note: not normalized!
inline double4 quat_axis(quat q) {
	return {q.x, q.y, q.z, 0};
}

// TODO is there a loss of precision due to acos?
inline double quat_angle(quat q) {
	//s = q.x * q.x + q.y * q.y + q.z * q.z;
	//return 2 * atan2(sqrt(s), q.w);

	return acos(q.w) * 2;
}

// single both are assumed to be unit quats, this function returns
// cos half angle to rotate from a to b
inline double quat_dot(quat a, quat b) {
	return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

// first row of matrix
inline double4 quat_dir_x(quat q) {
	double x = q.x, y = q.y, z = q.z, w = q.w;
	return double4{0.5 - (y * y + z * z), x * y - w * z, x * z + w * y, 0} * 2;
}

// second row of matrix
inline double4 quat_dir_y(quat q) {
	double x = q.x, y = q.y, z = q.z, w = q.w;
	return double4{x * y + w * z, 0.5 - (x * x + z * z), y * z - w * x, 0} * 2;
}

// third row of matrix
inline double4 quat_dir_z(quat q) {
	double x = q.x, y = q.y, z = q.z, w = q.w;
	return double4{x * z - w * y, y * z + w * x, 0.5 - (x * x + y * y), 0} * 2;
}


// first col of matrix (or first row of inv quat, assuming unit quat)
inline double4 quat_idir_x(quat q) {
	double x = q.x, y = q.y, z = q.z, w = q.w;
	return double4{0.5 - (y * y + z * z), x * y + w * z, x * z - w * y, 0} * 2;
}

/*
	Vector3 idirY() {
	double x = q.x, y = q.y, z = q.z, w = q.w;
		return new Vector3(2 * (x * y - w * z), 1 - 2 * (x * x + z * z), 2 * (y * z + w * x));
	}

	Vector3 idirZ() {
	double x = q.x, y = q.y, z = q.z, w = q.w;
		return new Vector3(2 * (x * z + w * y), 2 * (y * z - w * x), 1 - 2 * (x * x + y * y));
	}*/

// rotate vector by quaternion, returns q * quat(a, 0) * inv(q)
inline double3 quat_rotate(quat q, double3 a) {
	double iw = q.x * a.x + q.y * a.y + q.z * a.z;
	double ix = q.w * a.x + q.y * a.z - q.z * a.y;
	double iy = q.w * a.y - q.x * a.z + q.z * a.x;
	double iz = q.w * a.z + q.x * a.y - q.y * a.x;

	double x = iw * q.x + ix * q.w - iy * q.z + iz * q.y;
	double y = iw * q.y + ix * q.z + iy * q.w - iz * q.x;
	double z = iw * q.z - ix * q.y + iy * q.x + iz * q.w;
	return {x, y, z};
}

// column matrix
inline mat33 quat_to_matrix(quat q) {
	double x = q.x * 2;
	double y = q.y * 2;
	double z = q.z * 2;
	double xx = q.x * x;
	double yy = q.y * y;
	double zz = q.z * z;
	double xy = q.x * y;
	double xz = q.x * z;
	double yz = q.y * z;
	double xw = q.w * x;
	double yw = q.w * y;
	double zw = q.w * z;

	double3 a = {1 - (yy + zz), xy + zw, xz - yw};
	double3 b = {xy - zw, 1 - (xx + zz), yz + xw};
	double3 c = {xz + yw, yz - xw, 1 - (xx + yy)};
	return {a, b, c};
}

// rotate vector by inverse quaternion, returns inv(q) * quat(a, 0) * q
inline double4 quat_irotate(quat q, double4 a) {
	double iw = q.x * a.x + q.y * a.y + q.z * a.z;
	double ix = q.w * a.x - q.y * a.z + q.z * a.y;
	double iy = q.w * a.y + q.x * a.z - q.z * a.x;
	double iz = q.w * a.z - q.x * a.y + q.y * a.x;

	double x = iw * q.x + ix * q.w + iy * q.z - iz * q.y;
	double y = iw * q.y - ix * q.z + iy * q.w + iz * q.x;
	double z = iw * q.z + ix * q.y - iy * q.x + iz * q.w;
	return {x, y, z};
}

// Linear Interpolation
inline quat lerp(quat p, quat q, double t) {
	return p * (1 - t) + q * t;
}

// Spherical Linear Interpolation
inline quat slerp(quat p, quat q, double t) {
	double d = dot(p, q);
	if (d >= 0) {
		double a = acos(d);
		double k = 1 / sin(a);
		return p * (sin(a - t * a) * k) + q * (sin(t * a) * k);
	}
	double a = acos(-d);
	double k = 1 / sin(a);
	return p * (sin(a - t * a) * k) - q * (sin(t * a) * k);
}

inline quat quat_from_rot_mat(mat33 m) {
	double m00 = m.a.x;
	double m11 = m.b.y;
	double m22 = m.c.z;

	double x = 1 + m00 - m11 - m22;
	double y = 1 - m00 + m11 - m22;
	double z = 1 - m00 - m11 + m22;
	double w = 1 + m00 + m11 + m22;
	quat q = {x, y, z, w};
	q = sqrt(vmax(0, q)) / 2;

	if (m.c.y < m.b.z)
		q.x = -q.x;
	if (m.a.z < m.c.x)
		q.y = -q.y;
	if (m.b.x < m.a.y)
		q.z = -q.z;
	return q;
}
