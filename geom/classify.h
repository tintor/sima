#pragma once
#include <geom/aabb.h>
#include <geom/mesh.h>
#include <geom/plane.h>
#include <geom/polygon.h>
#include <geom/segment.h>

inline aabb2 Box(const polygon2& p) { return aabb2(p).buffer(Tolerance); }
inline aabb2 Box(const xpolygon2& p) { return aabb2(p.vertices()).buffer(Tolerance); }
inline aabb3 Box(cspan<double3> s) { return aabb3(s).buffer(Tolerance); }

// All Sign() functions return:
// +1 above / positive side
//  0 contact / between
// -1 below / negative side

inline constexpr int Sign(double d, double eps = Tolerance) {
    if (d > eps) return 1;
    if (d < -eps) return -1;
    return 0;
}

// TESTED
inline int Sign(segment2 s, double2 v) { return Sign(signed_double_area(s.a, s.b, v) / length(s.a - s.b)); }

inline int Sign(ray2 s, double2 v) { return Sign(signed_double_area(s.origin, s.origin + s.unit_dir, v)); }

inline int Sign(plane p, double3 v) { return Sign(p.distance(v)); }

struct IContact2 {
    double2 normal;
    double2 sa, sb;
};

struct IContact {
    enum class Type { Point, Segment, Polygon };
    Type type;
    double3 normal;  // from mb to ma
    double3 sa, sb;
};

struct Feature2 {
    double2 sa, sb;
};

struct XContact2 {
    double2 common;
    Feature2 fa, fb;
};

// All Classify() functions return:
// +1 disjoint / separate
//  0 touching
// -1 overlap / penetrating

// 2D
// ==
inline int Classify(segment2 s, double2 v) { return (squared(s.nearest(v) - v) > squared(Tolerance)) ? 1 : 0; }

// TESTED polygon vs point
bool PointInPolygon(double2 p, cspan<double2> polygon);
int Classify(const polygon2& a, double2 p, aabb2 box);
int Classify(const xpolygon2& a, double2 p, aabb2 box);
inline int Classify(const polygon2& a, double2 p) { return Classify(a, p, Box(a)); }
inline int Classify(const xpolygon2& a, double2 p) { return Classify(a, p, Box(a)); }

// polygon vs segment
using dpair = pair<double, double>;
int Classify(const polygon2& f, segment2 s, aabb2 box, vector<dpair>* intersections = nullptr);
int Classify(const xpolygon2& f, segment2 s, aabb2 box, vector<dpair>* intersections = nullptr);
inline int Classify(const polygon2& f, segment2 s, vector<dpair>* intersections = nullptr) {
    return Classify(f, s, Box(f), intersections);
}
inline int Classify(const xpolygon2& f, segment2 s, vector<dpair>* intersections = nullptr) {
    return Classify(f, s, Box(f), intersections);
}

// polygon vs polygon
int Classify(const xpolygon2& a, const xpolygon2& b, vector<IContact2>* contacts = nullptr);
int Classify(const polygon2& a, const polygon2& b, vector<IContact2>* contacts = nullptr);

// polygon vs ray
int Classify(const xpolygon2& f, ray2 s);

// 3D
// ==

int Classify(const face& f, double3 v, const aabb3& box);
int Classify(const face& f, const ray3& s, double* travel = nullptr);
pair<int, int> ClassifyDoubleSided(const face& f, const ray3& s, const aabb3& box);
int Classify(plane p, const segment3& s);
int Classify(const face& f, const segment3& s, bool reverse, vector<dpair>* intersections = nullptr,
             vector<IContact>* contacts = nullptr);
int Classify(const xmesh3& m, double3 p, const aabb3& box);
int Classify(const xmesh3& m, const segment3& s, const aabb3& box, bool reverse, vector<IContact>* contacts = nullptr);
int Classify(const xmesh3& ma, const xmesh3& mb, vector<IContact>* contacts = nullptr);
