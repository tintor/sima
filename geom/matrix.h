#pragma once
#include <geom/vector.h>

// each field is col of matrix
struct double22 {
    double2 a, b;

	static constexpr double22 make(double e) { return {{e, e}, {e, e}}; }

	double22 operator*(double k) const { return { a * k, b * k }; }
	double2 operator*(double2 v) const { return a * v.x + b * v.y; }
	double22 operator*(double22 q) const { return { *this * q.a, *this * q.b }; }

	double22 operator+(double22 m) const { return {a + m.a, b + m.b}; }
	double22 operator-(double22 m) const { return {a - m.a, b - m.b}; }
	void operator+=(double22 m) { a += m.a; b += m.b; }
	void operator-=(double22 m) { a -= m.a; b -= m.b; }
};

constexpr double22 double22_zero = {{0, 0}, {0, 0}};
constexpr double22 double22_identity = {{1, 0}, {0, 1}};

// each field is col of matrix
struct double33 {
    double3 a, b, c;

	static constexpr double33 make(double e) { return {{e, e, e}, {e, e, e}, {e, e, e}}; }

	double33 operator*(double k) const { return { a * k, b * k, c * k }; }
	double3 operator*(double3 v) const { return a * v.x + b * v.y + c * v.z; }
	double33 operator*(double33 q) const { return { *this * q.a, *this * q.b, *this * q.c }; }

	double33 operator+(double33 m) const { return {a + m.a, b + m.b, c + m.c}; }
	double33 operator-(double33 m) const { return {a - m.a, b - m.b, c - m.c}; }
	void operator+=(double33 m) { a += m.a; b += m.b; c += m.c; }
	void operator-=(double33 m) { a -= m.a; b -= m.b; c -= m.c; }
};

constexpr double33 double33_zero = {{0, 0, 0}, {0, 0, 0}};
constexpr double33 double33_identity = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

// each field is col of matrix
struct double44 {
    double4 a, b, c, d;

	static constexpr double44 make(double e) { return {{e, e, e, e}, {e, e, e, e}, {e, e, e, e}, {e, e, e, e}}; }

	double44 operator*(double k) const { return { a * k, b * k, c * k, d * k }; }
	double4 operator*(double4 v) const { return a * v.x + b * v.y + c * v.z + d * v.w; }
	double44 operator*(double44 q) const { return { *this * q.a, *this * q.b, *this * q.c, *this * q.d }; }

	double44 operator+(double44 m) const { return {a + m.a, b + m.b, c + m.c, d + m.d}; }
	double44 operator-(double44 m) const { return {a - m.a, b - m.b, c - m.c, d - m.d}; }
	void operator+=(double44 m) { a += m.a; b += m.b; c += m.c; d += m.d; }
	void operator-=(double44 m) { a -= m.a; b -= m.b; c -= m.c; d -= m.d; }
};

constexpr double44 double44_zero = {{0, 0, 0, 0}, {0, 0, 0, 0}};
constexpr double44 double44_identity = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

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

inline double trace(double22 m) { return m.a.x + m.b.y; }
inline double trace(double33 m) { return m.a.x + m.b.y + m.c.z; }
inline double trace(double44 m) { return m.a.x + m.b.y + m.c.z + m.d.w; }

inline double det(double22 m) { return det(m.a, m.b); }
inline double det(double33 m) { return det(m.a, m.b, m.c); }
inline double det(double44 m) { return det(m.a, m.b, m.c, m.d); }
