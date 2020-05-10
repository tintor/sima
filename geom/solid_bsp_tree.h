#pragma once
#include <geom/mesh.h>
#include <geom/plane.h>
#include <geom/triangle.h>

#include <random>

class SolidBSPTree {
   public:
    SolidBSPTree(const mesh3& mesh, uint num_samples, std::default_random_engine& rnd);
    bool intersects(double3 v);

   private:
    struct Node {
        plane divider;
        uint positive, negative;
    };

    struct Hist {
        int positive;
        int negative;
        int overlap;
        int stradle;
    };

    struct itriangle {
        ushort a, b, c;
    };

    struct BuildData {
        vector<double3> vertices;
        vector<itriangle> ifaces;
        vector<triangle3> faces;
        vector<double3> samples;
        std::default_random_engine rnd;
    };

    pair<uint, float> build_internal(BuildData& data, const vector<uint>& mesh, uint* samples_begin, uint* samples_end);

    void print_tree(uint n, int depth = 0);

    static void evaluate_candidate(plane candidate, BuildData& data, const vector<uint>& mesh, uint* samples_begin,
                                   uint* samples_end, long& best_heuristic, Hist& best_hist, plane& best_candidate);

    uint add_node() {
        uint n = node.size();
        node.resize(n + 1);
        percent.resize(n + 1);
        hist.resize(n + 1);
        box_size.resize(n + 1);
        return n;
    }

    vector<Node> node;
    vector<float> percent;
    vector<Hist> hist;
    vector<double3> box_size;
};
