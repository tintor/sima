#pragma once
#include "triangle.h"
#include "ivec.h"
#include <random>

struct iplane3 {
	ivec3 normal;
	int d;

	iplane3() { }

	iplane3(ivec3 normal, int d) : normal(normal), d(d) { }

	iplane3(const itriangle3& m) {
		THROW(not_implemented);
	}

	long distance(ivec3 v) const {
		return addi(doti(v, normal), d);
	}

	static int sign(ivec3 a, ivec3 b, ivec3 c, ivec3 v) {
		lvec3 normal = normali(a, b, c);
		return ::sign(doti(vconvert(normal, llvec3), vconvert(subi(v, a), llvec3)));
	}
};

class SolidBSPTree {
public:
    SolidBSPTree(const imesh3& mesh, uint num_samples, std::default_random_engine& rnd);
    bool intersects(ivec3 v);

private:
    struct Node {
        iplane3 divider;
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
        std::vector<ivec3> vertices;
        std::vector<itriangle> ifaces;
        std::vector<itriangle3> faces;
        std::vector<ivec3> samples;
        std::default_random_engine rnd;
    };

    std::pair<uint, float> build_internal(
            BuildData& data, const std::vector<uint>& mesh,
            uint* samples_begin, uint* samples_end);

    void print_tree(uint n, int depth = 0);

    static void evaluate_candidate(
            const iplane3& candidate,
            BuildData& data,
            const std::vector<uint>& mesh,
            uint* samples_begin,
            uint* samples_end,
            long& best_heuristic,
            Hist& best_hist,
            iplane3& best_candidate);

    uint add_node() {
        uint n = node.size();
        node.resize(n + 1);
        percent.resize(n + 1);
        hist.resize(n + 1);
        box_size.resize(n + 1);
        return n;
    }

    std::vector<Node> node;
    std::vector<float> percent;
    std::vector<Hist> hist;
    std::vector<ivec3> box_size;
};
