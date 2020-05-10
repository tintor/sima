#include <core/timestamp.h>
#include <geom/convex_hull.h>
#include <geom/is_valid.h>
#include <geom/mesh_import.h>
#include <geom/properties.h>
#include <geom/solid_bsp_tree.h>

#include <catch.hpp>

TEST_CASE("solid_bsp_tree", "[!hide][solid_bsp_tree]") {
    try {
        mesh3 mm = load_stl("models/bunny.stl");
        vector<double3> vertices;
        for (const triangle3& f : mm)
            for (auto i : range(3)) vertices.push_back(f[i]);
        print("hull...\n");
        mesh3 ch = convex_hull(vertices);
        print("IsValid %s\n", static_cast<int>(IsValid(mm)));
        print("Volume %s\n", Volume(mm));
        print("CenterOfMass %s\n", CenterOfMass(mm));
        print("IsConvex %s\n", is_convex(mm));

        print("IsValid %s\n", static_cast<int>(IsValid(ch)));
        print("Volume %s\n", Volume(ch));
        print("CenterOfMass %s\n", CenterOfMass(ch));
        print("IsConvex %s\n", is_convex(ch));

        std::default_random_engine rnd(0);

        Timestamp ta;
        SolidBSPTree tree(mm, 100000, rnd);
        print("%s\n", ta.elapsed_ms());
    } catch (std::exception& e) {
        print("std::runtime_error %s\n", e.what());
    }
}
