#include <core/exception.h>
#include <geom/aabb.h>
#include <geom/project.h>

ray2 Project(const ray3& s, int axis) { THROW(not_implemented); }

segment2 Project(const segment3& s, int axis) {
    if (axis == 0) return segment2(double2{s.a.y, s.a.z}, double2{s.b.y, s.b.z});
    if (axis == 1) return segment2(double2{s.a.x, s.a.z}, double2{s.b.x, s.b.z});
    if (axis == 2) return segment2(double2{s.a.x, s.a.y}, double2{s.b.x, s.b.y});
    THROW(invalid_argument);
}

double2 Project(double3 v, int axis) {
    if (axis == 0) return double2{v.y, v.z};
    if (axis == 1) return double2{v.x, v.z};
    if (axis == 2) return double2{v.x, v.y};
    THROW(invalid_argument);
}

// TODO issue this projection will stretch some axis more than others and skew tolerances!
xpolygon2 Project(const face& f, int axis) {
    size_t s = 0;
    for (auto ring : f) s += ring.size();
    xpolygon2 poly;
    poly.reserve(f.size(), s);

    if (axis == 0)
        for (auto ring : f) {
            poly.add({});
            for (auto p : ring) poly.add(double2{p.y, p.z});
        }
    if (axis == 1)
        for (auto ring : f) {
            poly.add({});
            for (auto p : ring) poly.add(double2{p.x, p.z});
        }
    if (axis == 2)
        for (auto ring : f) {
            poly.add({});
            for (auto p : ring) poly.add(double2{p.x, p.y});
        }
    return poly;
}

int ProjectionAxis(const face& f) {
    double3 s = aabb3(f.vertices()).size();
    if (s.x <= s.y && s.x <= s.z) return 0;
    if (s.y <= s.x && s.y <= s.z) return 1;
    return 2;
}
