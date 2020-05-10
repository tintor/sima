#pragma once
#include <core/format.h>
#include <core/span.h>
#include <geom/vector.h>

class sphere {
   public:
    sphere(double3 center, double radius) { s = extend(center, radius); }
    double3 center() const { return s.xyz; }
    double radius() const { return s.w; }
    bool contains(double3 p) const { return squared(center() - p) <= squared(radius()); }

   private:
    double4 s;
};

sphere minimal_sphere(sphere a, sphere b);
sphere minimal_sphere(sphere a, double3 b);
sphere minimal_sphere(double3 a, double3 b);
sphere minimal_sphere(double3 a, double3 b, double3 c);
sphere minimal_sphere(double3 a, double3 b, double3 c, double3 d);

sphere minimal_sphere(cspan<double3> points);

// not minimal, but close to it and faster to compute than minimal
sphere bounding_sphere(cspan<double3> points);
