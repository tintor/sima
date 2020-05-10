#pragma once

#include <geom/convex_hull.h>
#include <geom/mesh.h>
#include <geom/polygon.h>
#include <geom/primitives.h>

mesh3 generate_box(double3 size);

inline mesh3 generate_box(double sx, double sy, double sz) { return generate_box(double3{sx, sy, sz}); }

template <typename RND>
mesh3 generate_sphere(uint vertices, double radius, RND& rnd) {
    // generate N random vertices on sphere
    vector<double3> V;
    V.resize(vertices);
    for (auto i : range(vertices)) V[i] = uniform_dir3(rnd);

    // increase the distance between the closest two vertices
    // (as random clumps vertices together)
    vector<double3> delta;
    delta.resize(vertices);
    for (auto e : range(40)) {
        for (double3& v : delta) v = {0, 0, 0};
        for (auto i : range(vertices))
            for (auto j : range(i)) {
                double3 d = V[i] - V[j];
                d *= 0.05 / dot(d, d);
                delta[i] += d;
                delta[j] -= d;
            }
        for (auto i : range(vertices)) V[i] = normalize(V[i] + delta[i]);
    }

    vector<double3> I;
    I.resize(vertices);
    for (auto i : range(vertices)) I[i] = V[i] * radius;

    return convex_hull(I);
}

mesh3 generate_cross(double inner, double outer);

// cone (rmin or rmax = 0) or cylinder
mesh3 generate_cylinder(uint sides, double rmin, double rmax, double zmin, double zmax);

// platonic solids
mesh3 generate_tetrahedron(double radius);
mesh3 generate_cube(double radius);
mesh3 generate_octahedron(double radius);
mesh3 generate_dodecahedron(double radius);
mesh3 generate_icosahedron(double radius);

mesh3 generate_prism(const polygon2& poly, double zmin, double zmax);

mesh3 generate_pipe(cspan<double3> path, double radius, double vertices);

mesh3 generate_regular_polyhedra2(cspan<pair<int, int>> faces);
mesh3 generate_regular_polyhedra(cspan<cspan<int>> faces);
