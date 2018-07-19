#ifndef __SHAPE_HH__
#define __SHAPE_HH__

#include <vector>
#include <string>
#include <random>
#include "rendering.hh"
#include "primitives.hh"

typedef std::vector<triangle3> Mesh3d;

Mesh3d load_stl(const std::string& filename);
Mesh3d load_ply(const std::string& filename);

class SolidBSPTree {
public:
    SolidBSPTree(const Mesh3d& mesh, uint32_t num_samples, std::default_random_engine& rnd);
    bool intersects(const dvec3& v);

private:
    struct Node {
        plane divider;
        uint32_t positive, negative;
    };

    struct Hist {
        int positive;
        int negative;
        int overlap;
        int stradle;
    };

    struct itriangle {
        uint16_t a, b, c;
    };

    struct BuildData {
        std::vector<dvec3> vertices;
        std::vector<itriangle> ifaces;
        std::vector<triangle3> faces;
        std::vector<dvec3> samples;
        std::default_random_engine rnd;
    };

    std::pair<uint32, real> build_internal(
            BuildData& data, const std::vector<uint32_t>& mesh,
            uint32_t* samples_begin, uint32_t* samples_end);

    void print(uint32_t n, int depth = 0);

    static void evaluate_candidate(
            const plane& candidate,
            BuildData& data,
            const std::vector<uint32_t>& mesh,
            uint32_t* samples_begin,
            uint32_t* samples_end,
            int64_t& best_heuristic,
            Hist& best_hist,
            plane& best_candidate);

    uint32_t add_node() {
        uint32_t n = node.size();
        node.resize(n + 1);
        percent.resize(n + 1);
        hist.resize(n + 1);
        box_size.resize(n + 1);
        return n;
    }

    std::vector<Node> node;
    std::vector<float> percent;
    std::vector<Hist> hist;
    std::vector<dvec3> box_size;
};

// Mesh properties
// ===============

enum class Validity {
	OK = 0,
	TooFewFaces = 1,
	EdgeTooShort = 2,
	OpenEdge = 3,
	SeparateComponents = 4,
	SelfIntersection = 5,
};

Validity is_valid(const Mesh3d& mesh);

// Volume of a valid polyhedron
real volume(const Mesh3d& mesh);

// Center of mass of a valid polyhedron
dvec3 center_of_mass(const Mesh3d& mesh);

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
dmat3 moment_of_inertia(const Mesh3d& mesh);

// Is this valid mesh a convex polyhedron?
bool is_convex(const Mesh3d& mesh);

// Mesh generators
// ===============

// Returns empty vector if no solution (points are coplanar)
Mesh3d build_convex_hull(const std::vector<dvec3>& points);

Mesh3d BoxMesh(real sx, real sy, real sz);

template<typename RND>
Mesh3d SphereMesh(int vertices, RND& rnd) {
	std::vector<dvec3> V(vertices);
	FOR(i, vertices)
		V[i] = random_direction(rnd);
	return build_convex_hull(V);
}

Mesh3d CreateCrossMesh(real inner, real outer);

// Shape
// =====

// Shape is 3d solid, immutable, purely geometric and with origin in center of mass
class Shape {
public:
	Shape(const Mesh3d& mesh, /*in/out*/transform3& position);

	const auto& faces() const { return m_mesh; }
	const auto& convex_edges() const { return m_convex_edges; }
	const auto& convex_vertices() const { return m_convex_vertices; }
	real sphere_radius() const { return m_sphere_radius; }
	bool is_convex() const { return m_is_convex; }
	bool is_box() const { return m_is_box; }
	std::pair<dvec3, dvec3> box() const { return m_box; }

private:
	Mesh3d m_mesh;
	std::vector<segment3> m_convex_edges;
	std::vector<dvec3> m_convex_vertices;

	dmat3 m_inertia_tensor; // assuming density of 1kg/m^3
	real m_volume;
	real m_sphere_radius;
	std::pair<dvec3, dvec3> m_box; // oriented bounding box
	bool m_is_convex;
	bool m_is_box;
};

real signed_distance(const dvec3& v, const Shape& shape);

// either O(n*logn) SolidLeafBSPTree or O(n) algorithm using nearest face
inline bool intersects(const dvec3& v, const Shape& shape) {
	// TODO bounding box / sphere check
	// TODO optimize for box
	return signed_distance(v, shape) < 0;
}

real squared_distance_segment_origin(const segment3& p);

// TODO find all faces within epsilon distance from the edge
//      if any triangle is penetrated internally return true
//      if
//		rare_case:
//      collect at most one intersection point of edge with every nearby triangle
//		sort these points and check every mid point against the interior of shape
bool intersects_edge_interior(const segment3& edge, const Shape& shape);

bool intersects(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b);

struct Contact {
	real squared_dist;
	dvec3 position; // mid point
	dvec3 normal; // unit normal from body B to body A
};

bool is_vertex_triangle_contact(const dvec3& p, const triangle3& m, Contact& /*out*/contact);
bool is_edge_edge_contact(const segment3& p, const segment3& q, Contact& /*out*/contact);
// Assuming shapes aren't intersecting, look for pairs of features that are within ContactEpsilon
std::vector<Contact> find_all_contacts(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b);

// returns 0 if interecting
real distance(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b);

#endif