#pragma once
#include <geom/vector.h>

// each field is col of matrix
struct double22 {
    double2 a, b;
};

// each field is col of matrix
struct double33 {
    double3 a, b, c;
};

// each field is col of matrix
struct double44 {
    double4 a, b, c, d;
};

constexpr double22 IDENTITY_MAT22 = double22{ {1, 0}, {0, 1} };

inline double2 mul(double22 m, double2 v) {
	return m.a * v.x + m.b * v.y;
}

inline double3 mul(double33 m, double3 v) {
	return m.a * v.x + m.b * v.y + m.c * v.z;
}

inline double4 mul(double44 m, double4 v) {
	return m.a * v.x + m.b * v.y + m.c * v.z + m.d * v.w;
}

inline double22 transpose(double22 m) {
	return {{ m.a.x, m.b.x },
	        { m.a.y, m.b.y }};
}

inline double33 transpose(double33 m) {
	return {{ m.a.x, m.b.x, m.c.x },
	        { m.a.y, m.b.y, m.c.y },
	        { m.a.z, m.b.z, m.c.z }};
}

inline double44 transpose(double44 m) {
	return {{ m.a.x, m.b.x, m.c.x, m.d.x },
	        { m.a.y, m.b.y, m.c.y, m.d.y },
	        { m.a.z, m.b.z, m.c.z, m.d.z },
	        { m.a.w, m.b.w, m.c.w, m.d.w }};
}

inline double22 inv(double22 m) {
	double22 res;
    res.a.x = m.b.y;
    res.b.x = -m.b.x;
    res.a.y = -m.a.y;
    res.b.y = m.a.x;

	double det = 1 / (m.a.x * res.a.x + m.a.y * res.b.x);
	res.a *= det;
	res.b *= det;
    return res;
}

inline double33 inv(double33 m) {
	double33 res;
    res.a.x = m.b.y * m.c.z - m.c.y * m.b.z;
    res.b.x = -m.b.x * m.c.z + m.c.x * m.b.z;
    res.c.x = m.b.x * m.c.y - m.c.x * m.b.y;
    res.a.y = -m.a.y * m.c.z + m.c.y * m.a.z;
    res.b.y = m.a.x * m.c.z - m.c.x * m.a.z;
    res.c.y = -m.a.x * m.c.y + m.c.x * m.a.y;
    res.a.z = m.a.y * m.b.z - m.b.y * m.a.z;
    res.b.z = -m.a.x * m.b.z + m.b.x * m.a.z;
    res.c.z = m.a.x * m.b.y - m.b.x * m.a.y;

	double det = 1 / (m.a.x * res.a.x + m.a.y * res.b.x + m.a.z * res.c.x);
	res.a *= det;
	res.b *= det;
	res.c *= det;
    return res;
}

inline double44 inv(double44 m) {
	double44 res;

    res.a.x = m.b.y * m.c.z * m.d.w -
             m.b.y * m.c.w * m.d.z -
             m.c.y * m.b.z * m.d.w +
             m.c.y * m.b.w * m.d.z +
             m.d.y * m.b.z * m.c.w -
             m.d.y * m.b.w * m.c.z;

    res.b.x = -m.b.x * m.c.z * m.d.w +
              m.b.x * m.c.w * m.d.z +
              m.c.x * m.b.z * m.d.w -
              m.c.x * m.b.w * m.d.z -
              m.d.x * m.b.z * m.c.w +
              m.d.x * m.b.w * m.c.z;

    res.c.x = m.b.x * m.c.y * m.d.w -
             m.b.x * m.c.w * m.d.y -
             m.c.x * m.b.y * m.d.w +
             m.c.x * m.b.w * m.d.y +
             m.d.x * m.b.y * m.c.w -
             m.d.x * m.b.w * m.c.y;

    res.d.x = -m.b.x * m.c.y * m.d.z +
               m.b.x * m.c.z * m.d.y +
               m.c.x * m.b.y * m.d.z -
               m.c.x * m.b.z * m.d.y -
               m.d.x * m.b.y * m.c.z +
               m.d.x * m.b.z * m.c.y;

    res.a.y = -m.a.y * m.c.z * m.d.w +
              m.a.y * m.c.w * m.d.z +
              m.c.y * m.a.z * m.d.w -
              m.c.y * m.a.w * m.d.z -
              m.d.y * m.a.z * m.c.w +
              m.d.y * m.a.w * m.c.z;

    res.b.y = m.a.x * m.c.z * m.d.w -
             m.a.x * m.c.w * m.d.z -
             m.c.x * m.a.z * m.d.w +
             m.c.x * m.a.w * m.d.z +
             m.d.x * m.a.z * m.c.w -
             m.d.x * m.a.w * m.c.z;

    res.c.y = -m.a.x * m.c.y * m.d.w +
              m.a.x * m.c.w * m.d.y +
              m.c.x * m.a.y * m.d.w -
              m.c.x * m.a.w * m.d.y -
              m.d.x * m.a.y * m.c.w +
              m.d.x * m.a.w * m.c.y;

    res.d.y = m.a.x * m.c.y * m.d.z -
              m.a.x * m.c.z * m.d.y -
              m.c.x * m.a.y * m.d.z +
              m.c.x * m.a.z * m.d.y +
              m.d.x * m.a.y * m.c.z -
              m.d.x * m.a.z * m.c.y;

    res.a.z = m.a.y * m.b.z * m.d.w -
             m.a.y * m.b.w * m.d.z -
             m.b.y * m.a.z * m.d.w +
             m.b.y * m.a.w * m.d.z +
             m.d.y * m.a.z * m.b.w -
             m.d.y * m.a.w * m.b.z;

    res.b.z = -m.a.x * m.b.z * m.d.w +
              m.a.x * m.b.w * m.d.z +
              m.b.x * m.a.z * m.d.w -
              m.b.x * m.a.w * m.d.z -
              m.d.x * m.a.z * m.b.w +
              m.d.x * m.a.w * m.b.z;

    res.c.z = m.a.x * m.b.y * m.d.w -
              m.a.x * m.b.w * m.d.y -
              m.b.x * m.a.y * m.d.w +
              m.b.x * m.a.w * m.d.y +
              m.d.x * m.a.y * m.b.w -
              m.d.x * m.a.w * m.b.y;

    res.d.z = -m.a.x * m.b.y * m.d.z +
               m.a.x * m.b.z * m.d.y +
               m.b.x * m.a.y * m.d.z -
               m.b.x * m.a.z * m.d.y -
               m.d.x * m.a.y * m.b.z +
               m.d.x * m.a.z * m.b.y;

    res.a.w = -m.a.y * m.b.z * m.c.w +
              m.a.y * m.b.w * m.c.z +
              m.b.y * m.a.z * m.c.w -
              m.b.y * m.a.w * m.c.z -
              m.c.y * m.a.z * m.b.w +
              m.c.y * m.a.w * m.b.z;

    res.b.w = m.a.x * m.b.z * m.c.w -
             m.a.x * m.b.w * m.c.z -
             m.b.x * m.a.z * m.c.w +
             m.b.x * m.a.w * m.c.z +
             m.c.x * m.a.z * m.b.w -
             m.c.x * m.a.w * m.b.z;

    res.c.w = -m.a.x * m.b.y * m.c.w +
               m.a.x * m.b.w * m.c.y +
               m.b.x * m.a.y * m.c.w -
               m.b.x * m.a.w * m.c.y -
               m.c.x * m.a.y * m.b.w +
               m.c.x * m.a.w * m.b.y;

    res.d.w = m.a.x * m.b.y * m.c.z -
              m.a.x * m.b.z * m.c.y -
              m.b.x * m.a.y * m.c.z +
              m.b.x * m.a.z * m.c.y +
              m.c.x * m.a.y * m.b.z -
              m.c.x * m.a.z * m.b.y;

	double det = 1 / (m.a.x * res.a.x + m.a.y * res.b.x + m.a.z * res.c.x + m.a.w * res.d.x);
	res.a *= det;
	res.b *= det;
	res.c *= det;
	res.d *= det;
    return res;
}

inline double22 mul(const double22& p, const double22& q) {
	double22 e;
	e.a = mul(p, q.a);
	e.b = mul(p, q.b);
	return e;
}

inline double33 mul(const double33& p, const double33& q) {
	double33 e;
	e.a = mul(p, q.a);
	e.b = mul(p, q.b);
	e.c = mul(p, q.c);
	return e;
}

inline double44 mul(const double44& p, const double44& q) {
	double44 e;
	e.a = mul(p, q.a);
	e.b = mul(p, q.b);
	e.c = mul(p, q.c);
	e.d = mul(p, q.d);
	return e;
}

inline double det(double22 m) { return det(m.a, m.b); }
inline double det(double33 m) { return det(m.a, m.b, m.c); }
inline double det(double44 m) { return det(m.a, m.b, m.c, m.d); }
