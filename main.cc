#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <functional>
#include <stdexcept>
#include <cstdint>
#include <atomic>
#include <unistd.h>
#include <execinfo.h>

#define GLFW_INCLUDE_GLCOREARB
#ifndef __APPLE_CC__
    #include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>
#include "tinyformat.h"
#include "auto.h"
#include "rendering.hh"

using std::pair;
using std::cerr;
using std::cout;
using std::endl;
using std::size_t;
using std::min;
using std::max;

// ============

#include "primitives.hh"

typedef std::vector<triangle3> Mesh3d;

// ============

inline int64_t rdtsc() {
        uint lo, hi;
        __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
        return (static_cast<uint64_t>(hi) << 32) | lo;
}

struct Timestamp {
    Timestamp() : m_ticks(rdtsc()) { }

    int64_t elapsed(Timestamp a = Timestamp()) const { return a.m_ticks - m_ticks; }
    double elapsed_ms(Timestamp a = Timestamp()) const { return elapsed(a) * milisec_per_tick; }

    static void init(dvec3 a, i64vec3 b, dvec3 c, i64vec3 d) {
        glm::dvec3 q = (c - a) / glm::dvec3(d - b);
        if (q.x > q.y) std::swap(q.x, q.y);
        if (q.y > q.z) std::swap(q.y, q.z);
        if (q.x > q.y) std::swap(q.x, q.y);
        milisec_per_tick = q.y * 1000;
    }

    static double milisec_per_tick;

private:
    int64_t m_ticks;
};

double Timestamp::milisec_per_tick = 0;


// ============

struct SolidLeafBSPTree {
	static const SolidLeafBSPTree* Inside;

	SolidLeafBSPTree() { }
	SolidLeafBSPTree(SolidLeafBSPTree&) = delete;
	SolidLeafBSPTree& operator=(SolidLeafBSPTree&) = delete;

	~SolidLeafBSPTree() {
		destroy(positive);
		destroy(negative);
	}

	static void destroy(SolidLeafBSPTree* p) {
		if (p != nullptr && p != Inside)
			delete p;
	}

	plane divider;
	// nullptr means outside space, 0x0001 means inside space
	SolidLeafBSPTree* positive = nullptr;
	SolidLeafBSPTree* negative = nullptr;
};

const SolidLeafBSPTree* SolidLeafBSPTree::Inside = reinterpret_cast<const SolidLeafBSPTree*>(1);

int classify(dvec3 v, const plane& p) {
	auto d = p.distance(v);
	if (d > PlanarEpsilon)
		return 1;
	if (d < -PlanarEpsilon)
		return -1;
	return 0;
}

struct SolidBSPTreeBuildData {
	std::vector<dvec3> vertices;
	std::vector<dvec3> samples;
};

// TODO add memoization / add persisted memoization
// Returns pair of tree node and average query depth across all samples
pair<SolidLeafBSPTree*, real> build_solid_leaf_bsp_tree_internal(const SolidBSPTreeBuildData& data, const std::vector<uint16_t>& mesh, std::vector<uint32_t>& samples) {
	if (mesh.empty())
		return pair<SolidLeafBSPTree*, real>(nullptr, 0.0);

	if (mesh.size() / 3 >= 100) {
		// use heuristics instead of exaustive
		// TODO dividing plane generation:
		//      1) from any 3 vertices
		//      2) from an existing face
		//      3) major axis aligned to some vertex
		// TODO dividing plane heuristic ranking:
		//      1) avoid straddle
		//      2) equal number of faces on both sides
		//      3) early out (one side empty)

	}

	std::unordered_set<uint16_t> vertices;
	FOR_EACH(i, mesh)
		vertices.insert(i);

	// All possible candidate dividers (includes optimal solution)
	std::vector<plane> candidates;
	FOR(a, vertices.size())
		FOR(b, a)
			FOR(c, b)
				candidates.push_back(plane(data.vertices[a], data.vertices[b], data.vertices[c]));
	vertices.clear();

	SolidLeafBSPTree* tree = new SolidLeafBSPTree;
	real best_score = 1e100;
	std::vector<uint16_t> pmesh, nmesh;
	std::vector<uint32_t> psamples, nsamples;
	FOR_EACH(candidate, candidates) {
		// classify faces against divider and build tree recursively
		// TODO if it can fit repartition mesh array in-place
		pmesh.clear();
		nmesh.clear();
		FOR(f, mesh.size() / 3) {
			int mn = 1, mx = -1;
			FOR(i, 3) {
				int d = classify(data.vertices[mesh[f * 3 + i]], candidate);
				mn = std::min(mn, d);
				mx = std::max(mx, d);
			}
			if (mx == 1)
				FOR(i, 3)
					pmesh.push_back(mesh[f * 3 + i]);
			if (mn == -1)
				FOR(i, 3)
					nmesh.push_back(mesh[f * 3 + i]);
		}
		// TODO repartition samples array in place instead
		psamples.clear();
		nsamples.clear();
		FOR_EACH(s, samples)
			(candidate.distance(data.samples[s]) > 0 ? psamples : nsamples).push_back(s);
		// TODO handle case when pmesh or nmesh is empty AND set them to EMPTY or FULL properly!
		auto presult = build_solid_leaf_bsp_tree_internal(data, pmesh, psamples);
		auto nresult = build_solid_leaf_bsp_tree_internal(data, nmesh, nsamples);
		real score = 1 + (psamples.size() * presult.second + nsamples.size() * nresult.second) / samples.size();

		if (score < best_score) {
			tree->divider = candidate;
			SolidLeafBSPTree::destroy(tree->positive);
			SolidLeafBSPTree::destroy(tree->negative);
			tree->positive = presult.first;
			tree->negative = nresult.first;
			best_score = score;
		} else {
			SolidLeafBSPTree::destroy(presult.first);
			SolidLeafBSPTree::destroy(nresult.first);
		}
	}
	return pair(tree, best_score);
}

std::unique_ptr<SolidLeafBSPTree> build_solid_leaf_bsp_tree(const Mesh3d& mesh, uint32_t num_samples) {
	SolidBSPTreeBuildData data;

	// Init samples
	data.samples.resize(num_samples);
	std::vector<uint32_t> isamples(num_samples);
	std::default_random_engine rnd;
	FOR(i, num_samples) {
		data.samples[i] = random_vector(rnd); // TODO make sure inside bounding box / sphere
		isamples[i] = i;
	}

	// Init vertices and imesh
	std::unordered_map<dvec3, uint16_t> vec_index;
	std::vector<uint16_t> imesh;
	FOR_EACH(face, mesh)
		FOR(i, 3) {
			const auto& v = face[i];
			if (vec_index.count(v) == 0) {
				vec_index[v] = data.vertices.size();
				data.vertices.push_back(v);
			}
			imesh.push_back(vec_index[v]);
		}

	auto result = build_solid_leaf_bsp_tree_internal(data, imesh, isamples);
	return std::unique_ptr<SolidLeafBSPTree>(result.first);
}

// Note: may return either true or false for boundary!
bool intersects(const SolidLeafBSPTree* tree, const dvec3& v) {
	while (true) {
		if (tree == nullptr)
			return false;
		if (tree == SolidLeafBSPTree::Inside)
			return true;
		tree = (tree->divider.distance(v) > 0) ? tree->positive : tree->negative;
	}
}

// ============

// Valid if triangles are not intersecting, except in one shared edge or one shared vertex
bool are_valid_mesh_faces(const triangle3& a, const triangle3& b) {
	// Axis check for early exit
	dvec3 amin = mini(a), amax = maxi(a);
	dvec3 bmin = mini(b), bmax = maxi(b);
	FOR(i, 3)
		if (DisjointIntervals(amin[i], amax[i], bmin[i], bmax[i]))
			return true;

	// Plane check
	if (!intersects(a, plane(b)) || !intersects(b, plane(a)))
		return true;

	// Non-planar case
	int match[3] = {-1, -1, -1};
	FOR(i, 3) FOR(j, 3)
		if (a[i] == b[j])
			match[i] = j;
	int count = 0;
	FOR(i, 3)
		if (match[i] >= 0)
			count += 1;
	if (count == 1) {

	}

	// Planar case
	// TODO there is some intersection:
	// OK case is if intersection is one vertex of both A and B (and no overlap)
	// OK case is if intersection is one edge of both A and B (and no overlap)
	return true;
}

enum class Validity {
	OK = 0,
	TooFewFaces = 1,
	EdgeTooShort = 2,
	OpenEdge = 3,
	SeparateComponents = 4,
	SelfIntersection = 5,
};

Validity is_valid(const Mesh3d& mesh) {
	if (mesh.size() < 4)
		return Validity::TooFewFaces;

	// Minimal length of any edge is 10xPlanarEpsilon
	FOR(i, mesh.size())
		FOR_EACH_EDGE(a, b, mesh[i])
			if (squared(*a - *b) <= squared(10 * PlanarEpsilon))
				return Validity::EdgeTooShort;

	// TODO any two different vertices can't be closer than PlanarEpsilon
	// TODO any vertex can't be closer than PlanarEpsilon to any other edge

	// Extract all edges
	std::unordered_map<segment3, int> edge_to_face;
	FOR(i, mesh.size())
		FOR_EACH_EDGE(p, q, mesh[i])
 			edge_to_face[segment3(*p, *q)] = i;

	// Every edge must appear exactly twice (in opposite directions)
	FOR_EACH(e, edge_to_face)
		if (edge_to_face.count(e.first.reverse()) == 0)
			return Validity::OpenEdge;

	// All triangles must form a single component
	std::vector<UnionFind> component(mesh.size());
	FOR(i, mesh.size())
		FOR_EACH_EDGE(p, q, mesh[i])
			component[edge_to_face[segment3(*q, *p)]].Union(component[i]);
	FOR(i, component.size())
		if (component[i].Find() != component[0].Find())
			return Validity::SeparateComponents;

	// No self intersections
	FOR(i, mesh.size())
		FOR(j, i)
			if (!are_valid_mesh_faces(mesh[i], mesh[j]))
				return Validity::SelfIntersection;

	return Validity::OK;
}

/*TEST_CASE("is_valid") {
	REQUIRE(is_valid(std::vector<triangle3>{}) == Validity::TooFewFaces);
	dvec3 a(0, 0, 0);
	dvec3 b(1, 0, 0);
	dvec3 c(0, 1, 0);
	dvec3 d(0, 0, 1);
}*/

// Mesh properties
// ===============

// Volume of a valid polyhedron
real volume(const Mesh3d& mesh) {
	real volume = 0;
	FOR_EACH(m, mesh) {
		real s = 0, z = 0;
		FOR_EACH_EDGE(a, b, m) {
			z += b->z;
			s += (a->y + b->y) * (a->x - b->x);
		}
		volume += z * s;
	}
	return volume / 6;
}

// TODO test with volume of cube (randomly rotated)

// Center of mass of a valid polyhedron
dvec3 center_of_mass(const Mesh3d& mesh) {
	dvec3 P(0, 0, 0);
	real V = 0;
	FOR_EACH(f, mesh) {
		real v = dot(f.a, cross(f.b, f.c));
		P += (f.a + f.b + f.c) * v;
		V += v;
	}
	return P / (V * 4);
}

// Moment of inertia of a valid polyhedron with Center of Mass = 0 and Density = 1
dmat3 moment_of_inertia(const Mesh3d& mesh) {
	constexpr real a = 1 / 60., b = 1 / 120.;
    constexpr real canon[] = { a, b, b, b, a, b, b, b, a};
	const dmat3 canonical = glm::make_mat3(canon);
    dmat3 C = dmat3(0); // covariance
	FOR_EACH(f, mesh) {
		dmat3 A;
		FOR(i, 3)
			glm::row(A, i) = f[i]; // TODO: or column(i) = ?
		C += (transpose(A) * canonical * A) * determinant(A);
	}
	return full_mat3(C[0][0] + C[1][1] + C[2][2]) - C; // C -> I
}

// Is this valid mesh a convex polyhedron?
bool is_convex(const Mesh3d& mesh) {
	// Extract all unique vertices
	std::unordered_set<dvec3> vertices;
	FOR_EACH(f, mesh)
		vertices.insert(f.begin(), f.begin() + 2);

	FOR_EACH(f, mesh) {
		plane p(f);
		FOR_EACH(v, vertices)
			if (p.distance(v) > PlanarEpsilon)
				return false;
	}
	return true;
}

// Mesh generators
// ===============

// Returns empty vector if no solution (points are coplanar)
Mesh3d build_convex_hull(const std::vector<dvec3>& points) {
	// First two points (A and B) on the hull (extremes on X axis)
	size_t ai = 0, bi = 0;
	FOR(i, points.size()) {
		const dvec3& p = points[i];
		if (p.x < points[ai].x)
			ai = i;
		if (p.x > points[bi].x)
			bi = i;
	}
	if (ai == bi)
		return Mesh3d();
	dvec3 a = points[ai], b = points[bi];

	// Third point C on the hull (furthest from line AB)
	segment3 line(a, b);
	real max_dist2 = 0;
	size_t ci = 0;
	FOR(i, points.size()) {
		const dvec3& p = points[i];
		real dist2 = line_point_squared_distance(a, b, p);
		if (dist2 > max_dist2) {
			ci = i;
			max_dist2 = dist2;
		}
	}
	if (max_dist2 < squared(PlanarEpsilon))
		return Mesh3d();
	dvec3 c = points[ci];

	// Fourth point D on the hull (furthest from plane ABC)
	plane plane(a, b, c);
	real max_dist = 0;
	size_t di = 0;
	FOR(i, points.size()) {
		const dvec3& p = points[i];
		real dist = plane.distance(p);
		if (abs(dist) > abs(max_dist)) {
			di = i;
			max_dist = dist;
		}
	}
	if (abs(max_dist) < PlanarEpsilon)
		return Mesh3d();
	dvec3 d = points[di];

	// Construct initial tetrahedron hull,
	// All faces are oriented facing outside with right hand rule.
	Mesh3d hull;
	hull.reserve(4);
	if (max_dist < 0) {
		hull.push_back(triangle3(a, b, c));
		hull.push_back(triangle3(b, a, d));
		hull.push_back(triangle3(c, b, d));
		hull.push_back(triangle3(a, c, d));
	} else {
		hull.push_back(triangle3(c, b, a));
		hull.push_back(triangle3(a, b, d));
		hull.push_back(triangle3(b, c, d));
		hull.push_back(triangle3(c, a, d));
	}

	// Expand hull to include all remaining points outside of it
	FOR(pi, points.size()) {
		const dvec3& p = points[pi];
		if (pi == ai || pi == bi || pi == ci || pi == di)
			continue;

		// Remove faces on hull covered by new vertex
		std::unordered_set<segment3> open_edges;
		FOR(i, hull.size()) {
			const triangle3& f = hull[i];
			dvec3 normal = Normal(f.a, f.b, f.c);
			real dist = dot(normal, p - f.a);

			// If P is in front of face F
			if (dist > 0 && dist * dist > squared(PlanarEpsilon) * squared(normal)) {
				// Add edges of removed face to open_edges
				FOR_EACH_EDGE(a, b, f) {
					segment3 e(*a, *b);
					// If two faces that share an edge are both removed,
					// then their shared edge isn't open anymore.
					if (open_edges.erase(e.reverse()) == 0)
						open_edges.insert(e);
				}

				// Remove face
				hull[i] = hull[hull.size() - 1];
				hull.resize(hull.size() - 1);
				i -= 1;
			}
		}

		// For each open edge create a face that connects it with P
		FOR_EACH(e, open_edges)
			hull.push_back(triangle3(e.a, e.b, p));
	}
	return hull;
}

Mesh3d BoxMesh(real sx, real sy, real sz) {
	std::vector<dvec3> vertices;
	vertices.reserve(8);
	for (int x = -1; x <= 1; x += 2)
		for (int y = -1; y <= 1; y += 2)
			for (int z = -1; z <= 1; z += 2)
				vertices.push_back(dvec3(x * sx, y * sy, z * sz));
	return build_convex_hull(vertices);
}

Mesh3d SphereMesh(int vertices) {
	std::default_random_engine rnd;
	std::vector<dvec3> V(vertices);
	FOR(i, vertices)
		V[i] = random_direction(rnd);
	return build_convex_hull(V);
}

/*TEST_CASE("SphereMesh - ConvexHull - IsConvex") {
	Mesh3d mesh = SphereMesh(100);
	real v = 4.0 / 3 * M_PI;
	REQUIRE(abs(volume(mesh) - v) / v < 0.15);
	REQUIRE(is_valid(mesh) == Validity::OK);
	REQUIRE(is_convex(mesh));
}*/

void AddQuad(Mesh3d& mesh, dvec3 a, dvec3 b, dvec3 c, dvec3 d) {
	mesh.push_back(triangle3(a, b, c));
	mesh.push_back(triangle3(c, d, a));
}

Mesh3d CreateCrossMesh(real inner, real outer) {
	Mesh3d m;
	// TODO finish
	real a = inner, b = outer;

	// square face - Z axis
	AddQuad(m, dvec3(a, a, b), dvec3(-a, a, b), dvec3(-a, -a, b), dvec3(a, -a, b));
	AddQuad(m, dvec3(a, a, -b), dvec3(-a, -a, -b), dvec3(-a, a, -b), dvec3(a, -a, -b));
	// square face - Y axis
	AddQuad(m, dvec3(a, b, a), dvec3(-a, b, a), dvec3(-a, b, -a), dvec3(a, b, -a));
	AddQuad(m, dvec3(a, -b, a), dvec3(-a, -a, -a), dvec3(-a, -b, a), dvec3(a, -b, -a));
	// square face - X axis
	AddQuad(m, dvec3(b, a, a), dvec3(b, -a, a), dvec3(b, -a, -a), dvec3(b, a, -a));
	AddQuad(m, dvec3(-b, a, a), dvec3(-b, -a, -a), dvec3(-b, -a, a), dvec3(-b, a, -a));
	return m;
}

/*TEST_CASE("CreateCrossMesh") {
	Mesh3d cross = CreateCrossMesh(0.05, 0.25);
	//REQUIRE(IsValid(cross) == Validity::OK);
}*/

bool lexicographical_less(dvec3 a, dvec3 b) {
	return a.x < b.x || (a.x == b.x && a.y < b.y) || (a.x == b.x && a.y == b.y && a.z < b.z);
}

template<typename Container, typename Func>
void sort(Container& container, const Func& func) {
	std::sort(container.begin(), container.end(), func);
}

// Shape is 3d solid, immutable, purely geometric and with origin in center of mass
class Shape {
public:
	Shape(const Mesh3d& mesh, /*in/out*/transform3& position) {
		if (is_valid(mesh) != Validity::OK)
			throw new std::runtime_error("invalid mesh");

		// Move to center of mass
		dvec3 center = center_of_mass(mesh);
		m_mesh = mesh; // copy
		FOR_EACH(face, m_mesh)
			face -= center;
		// TODO update 'position' out parameter, as shape is moved

		m_volume = volume(m_mesh);
		m_inertia_tensor = moment_of_inertia(m_mesh);

		// TODO HARD? rotate polyhedron, so inertia tensor only has primary axes (but return back its original orientation)
		//      so inertia tensor becomes inertia vector. Will this rotation make OBB smaller?

		m_is_convex = ::is_convex(m_mesh);

		// Compute convex edges and vertices
		std::unordered_map<segment3, dvec3> third_vertex;
		FOR_EACH(face, m_mesh)
			FOR(i, 3) {
				const dvec3& a = face[i], b = face[(i + 1) % 3], c = face[(i + 2) % 3];
				third_vertex[segment3(a, b)] = c;
			}
		std::unordered_set<dvec3> set_of_convex_vertices;
		FOR_EACH(face, m_mesh)
			FOR(i, 3) {
				const dvec3& a = face[i], b = face[(i + 1) % 3], c = face[(i + 2) % 3];
				if (edge_angle(a, b, c, third_vertex[segment3(b, a)]) < M_PI) {
					if (lexicographical_less(a, b))
						m_convex_edges.push_back(segment3(a, b));
				} else {
					set_of_convex_vertices.erase(a);
					set_of_convex_vertices.erase(b);
				}
			}
        FOR_EACH(v, set_of_convex_vertices)
			m_convex_vertices.push_back(v);

		// Compute radius of bounding sphere and size of object bounding box
		m_sphere_radius = 0;
		m_box.first = m_box.second = m_convex_vertices[0];
		FOR_EACH(v, m_convex_vertices) {
			m_sphere_radius = std::max(m_sphere_radius, squared(v));
			m_box.first = min(m_box.first, v);
			m_box.second = max(m_box.second, v);
		}
		m_sphere_radius = sqrt(m_sphere_radius);

		// Sort faces by decreasing surface area
		sort(m_mesh, [](const triangle3& p, const triangle3& q) {
			return p.squared_area_x4() > q.squared_area_x4();
		});
		// Sort convex edges by decreasing length
		sort(m_convex_edges, [](const segment3& p, const segment3& q) {
			return squared(p.a - p.b) > squared(q.a - q.b);
		});
		// Sort convex vertices by decreasing distance from center
		sort(m_convex_vertices, [](const dvec3& p, const dvec3& q) {
			return squared(p) > squared(q);
		});
	}

	const auto& faces() const { return m_mesh; }
	const auto& convex_edges() const { return m_convex_edges; }
	const auto& convex_vertices() const { return m_convex_vertices; }
	bool is_convex() const { return m_is_convex; }
	bool is_box() const { return false; }
	real sphere_radius() const { return m_sphere_radius; }
	pair<dvec3, dvec3> box() const { return m_box; }

private:
	Mesh3d m_mesh;
	std::vector<segment3> m_convex_edges;
	std::vector<dvec3> m_convex_vertices;

	dmat3 m_inertia_tensor; // assuming density of 1kg/m^3
	real m_volume;
	real m_sphere_radius;
	pair<dvec3, dvec3> m_box; // oriented bounding box
	bool m_is_convex;
};

/*const auto random_directions = [](){
	std::array<dvec3, 128> dirs;
	std::default_random_engine rnd;
	FOR_EACH(d, dirs)
		d = random_direction(rnd);
	return dirs;
}();*/

struct Body {
	Shape shape;
	dvec3 position;
	dquat orientation;
};

real signed_distance(const dvec3& v, const Shape& shape) {
	// Can't just take the sdist sign of min_dist as two faces that share the edge can have the same dist but different sgn(sdist)
	real min_dist = std::numeric_limits<real>::max(), max_sdist = 0;
	constexpr real eps = ContactEpsilon; // ??? arbitrary, need to account for rounding errors when computing distance(vertex, face)
	FOR_EACH(face, shape.faces()) {
		plane p(face); // TODO precompute this

		real dist, sdist;
		FOR_EACH_EDGE(a, b, face)
			if (plane::sign(*a, *b, *a + p.normal, v) > 0) {
				dist = distance(v, segment3(*a, *b));
				if (dist > min_dist + eps)
					goto next;

				sdist = p.distance(v);
				goto rest;
			}

		sdist = p.distance(v);
		dist = abs(sdist);
		if (dist > min_dist + eps)
			continue;
		rest:
		if (dist < min_dist) {
			if (dist + eps < min_dist)
				max_sdist = sdist;
			min_dist = dist;
		}
		if (abs(sdist) > abs(max_sdist))
			max_sdist = sdist;
		next:;
	}
	return max_sdist < 0 ? -min_dist : min_dist;
}

// either O(n*logn) SolidLeafBSPTree or O(n) algorithm using nearest face
bool intersects(const dvec3& v, const Shape& shape) {
	// TODO bounding box / sphere check
	// TODO optimize for box
	return signed_distance(v, shape) < 0;
}

real squared_distance_segment_origin(const segment3& p) {
	dvec3 d = p.b - p.a;
	real t = -dot(p.a, d);
	if (t <= 0)
		return squared(p.a);
	real dd = squared(d);
	if (t >= dd)
		return squared(p.b);
	return squared(p.a + d * (t / dd));
}

// TODO find all faces within epsilon distance from the edge
//      if any triangle is penetrated internally return true
//      if
//		rare_case:
//      collect at most one intersection point of edge with every nearby triangle
//		sort these points and check every mid point against the interior of shape
bool intersects_edge_interior(const segment3& edge, const Shape& shape) {
	// Quick bounding sphere check
	// TODO precompute right hand side and store in Shape
	if (squared_distance_segment_origin(edge) > squared(shape.sphere_radius()))
		return false;

	// Quick bounding box check
	//if (!is_edge_intersecting_box())
	//	return false;
	if (shape.is_box())
		return true;

	// Note: assumes edge vertices are outside of shape
	constexpr real eps = 1e-8; // TODO
	/*std::vector<triangle3> nearest;
	FOR_EACH(face, shape.faces()) {
		// if edge goes strictly through interior of triangle TODO and edge vertices are not on the plane
		if (intersects_in_point(edge, face)
				&& distance(edge, segment3(face[0], face[1])) > eps
				&& distance(edge, segment3(face[1], face[2])) > eps
				&& distance(edge, segment3(face[2], face[0])) > eps)
			return true;
		real dist = disjoint_distance(edge, face);
		if (dist < eps) {
			nearest.push_back(face);
		}
	}*/
	throw new std::runtime_error("unfinished");
}

pair<real, real> project_obb(dvec3 obb_position, dvec3 obb_size, dmat3 obb_orientation, dvec3 dir) {
    // TODO
	return pair<real, real>(0.0, 0.0);
}

bool are_oriented_boxes_intersecting(const pair<dvec3, dvec3>& box_a, const transform3& pos_a, const pair<dvec3, dvec3>& box_b, const transform3& pos_b) {
	throw new std::runtime_error("unfinished");
}

bool approx_intersects(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	// Bounding sphere check
	if (squared(pos_a.position - pos_b.position) > squared(shape_a.sphere_radius() + shape_b.sphere_radius()))
		return false;

	return are_oriented_boxes_intersecting(shape_a.box(), pos_a, shape_b.box(), pos_b);
}

bool intersects(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	// Check all vertices of one against the other
	FOR_EACH(vertex_a, shape_a.convex_vertices())
		if (intersects(pos_b.to_local(pos_a.to_global(vertex_a)), shape_b))
			return true;
	FOR_EACH(vertex_b, shape_b.convex_vertices())
		if (intersects(pos_a.to_local(pos_b.to_global(vertex_b)), shape_a))
			return true;

	// Check all edges of one against the other
	FOR_EACH(edge_a, shape_a.convex_edges())
		if (intersects_edge_interior(pos_b.to_local(pos_a.to_global(edge_a)), shape_b))
			return true;
	FOR_EACH(edge_b, shape_b.convex_edges())
		if (intersects_edge_interior(pos_a.to_local(pos_b.to_global(edge_b)), shape_a))
			return true;

	return false;
}

// TODO test intersects:
// - "simulate" two boxes "colliding" one another (different sizes and movement directions)
// - same test, but with bunnies
// - OpenGL visualization (for debugging the test)

struct Contact {
	real squared_dist;
	dvec3 position; // mid point
	dvec3 normal; // unit normal from body B to body A
};

bool is_vertex_triangle_contact(const dvec3& p, const triangle3& m, Contact& /*out*/contact) {
	dvec3 normal = Normal(m);

	// if P is outside of triangle
	FOR_EACH_EDGE(a, b, m)
		if (plane::sign(*a, *b, *a + normal, p) > 0)
			return false;

	// P is inside of triangle
	dvec3 nearest = p - normal * dot(normal, p - m.a);
	contact.squared_dist = squared(nearest - p);
	if (contact.squared_dist > squared(ContactEpsilon))
		return false;
	contact.position = (nearest + p) / static_cast<real>(2);
	contact.normal = glm::normalize(normal);
	return true;
}

bool is_edge_edge_contact(const segment3& p, const segment3& q, Contact& /*out*/contact) {
	dvec3 A = p.b - p.a, B = q.b - q.a, C = p.a - q.a;
	real aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);
	constexpr real tiny = 1e-8;

	real d = aa * bb - ab * ab;
	real s = ab * bc - bb * ac;
	real t = aa * bc - ab * ac;
	// Note: [d >= tiny * aa * bb] is needed to make parallelness test indepent of line lengths
	if ((d >= tiny && d >= tiny * aa * bb && 0 <= s && s <= d && 0 <= t && t <= d)
			|| (d <= -tiny && d <= -tiny * aa * bb && d <= s && s <= 0 && d <= t && t <= 0)) {
		dvec3 P = p.a + A * (s / d);
		dvec3 Q = q.a + B * (t / d);
		contact.squared_dist = squared(P - Q);
		if (contact.squared_dist > squared(ContactEpsilon))
			return false;
		contact.position = (P + Q) / static_cast<real>(2);
		// Note: using AxB here instead of P-Q as P and Q will be very close.
		contact.normal = glm::normalize(cross(A, B));
		if (dot(P - Q, contact.normal) < 0)
			contact.normal = -contact.normal;
		return true;
	}

	return false;
}

bool all_less_equal(const dvec3& a, const dvec3& b) {
	return a.x <= b.x && a.y <= b.y && a.z <= b.z;
}

bool approx_intersects(const dvec3& v, const Shape& shape) {
	return squared(v) <= squared(shape.sphere_radius()) && all_less_equal(shape.box().first, v) && all_less_equal(v, shape.box().second);
}

bool approx_intersects(const segment3& v, const Shape& shape) {
	// TODO
	return true;
}

// Assuming shares aren't intersecting, look for pairs of features that are within ContactEpsilon
std::vector<Contact> find_all_contacts(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	std::vector<Contact> contacts;
	Contact contact;

	// vertex vs face contacts
	FOR_EACH(vertex_a, shape_a.convex_vertices()) {
		dvec3 vertex_a_local = pos_b.to_local(pos_a.to_global(vertex_a));
		// TODO this approx check needs to be a little loose to allow for contacts
		if (!approx_intersects(vertex_a_local, shape_b))
			continue;
		FOR_EACH(face_b, shape_b.faces())
			if (is_vertex_triangle_contact(vertex_a_local, face_b, /*out*/contact)) {
				contact.position = pos_b.to_global(contact.position);
				contact.normal = pos_b.to_global_dir(contact.normal);
				contacts.push_back(contact);
			}
	}
	FOR_EACH(vertex_b, shape_b.convex_vertices()) {
		dvec3 vertex_b_local = pos_a.to_local(pos_b.to_global(vertex_b));
		// TODO this approx check needs to be a little loose to allow for contacts
		if (!approx_intersects(vertex_b_local, shape_a))
			continue;
		FOR_EACH(face_a, shape_a.faces())
			if (is_vertex_triangle_contact(vertex_b_local, face_a, /*out*/contact)) {
				contact.position = pos_a.to_global(contact.position);
				contact.normal = -pos_a.to_global_dir(contact.normal);
				contacts.push_back(contact);
			}
	}

	// edge vs edge contacts
	// TODO flip loop order if B has less edges than A
	FOR_EACH(edge_a, shape_a.convex_edges()) {
		segment3 edge_a_global = pos_a.to_global(edge_a);
		// TODO this approx check needs to be a little loose to allow for contacts
		if (!approx_intersects(pos_b.to_local(edge_a_global), shape_b))
			continue;
		// TODO avoid transforming every edge in B
		FOR_EACH(edge_b, shape_b.convex_edges())
			if (is_edge_edge_contact(edge_a_global, pos_b.to_global(edge_b), /*out*/contact))
				contacts.push_back(contact);
	}

	return contacts;
}

// returns 0 if interecting
real distance(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	// Upper bound on distance
    real dist = max(distance(pos_a.position, pos_b.position), shape_a.sphere_radius(), shape_b.sphere_radius());

	// vertex vs face
	FOR_EACH(vertex_a, shape_a.convex_vertices()) {
		dvec3 vertex_a_local = pos_b.to_local(pos_a.to_global(vertex_a));
		// TODO skip if current dist is smaller than distance between vertex_a and bounding sphere of B
		FOR_EACH(face_b, shape_b.faces())
            dist = std::min(dist, distance(face_b, vertex_a_local));
	}
	FOR_EACH(vertex_b, shape_b.convex_vertices()) {
		dvec3 vertex_b_local = pos_a.to_local(pos_b.to_global(vertex_b));
		// TODO skip if current dist is smaller than distance between vertex_b and bounding sphere of A
		FOR_EACH(face_a, shape_a.faces())
            dist = std::min(dist, distance(face_a, vertex_b_local));
	}

	// edge vs edge contacts
	FOR_EACH(edge_a, shape_a.convex_edges()) {
		segment3 edge_a_local = pos_b.to_local(pos_a.to_global(edge_a));
		// TODO skip if current dist is smaller than distance between edge_a and bounding sphere of B
		FOR_EACH(edge_b, shape_b.convex_edges())
            dist = std::min(dist, distance(edge_a_local, edge_b));
	}

    return dist;
}

// Dynamics and simulation
// =======================

struct DynamicBody : public Body {
	real mass;
	dvec3 velocity;
	dvec3 angular_velocity;
};

class Joint {
	// ball joint (common point) 3DF
	// hinge / axel joint (common edge) 1DF
	// cylindrical joint (common line) 2DF
	// prismatic joint 1DF (two common parallel lines)
};

class BallJoint {
	dvec3 pa, pb;
};

class HingeJoint {
	dvec3 pa, qa, pb, qb;
	dvec3 ra, rb; // just for computing relative angle (must be unit and normal to PQ)
};

// Open chain articulated only
// TODO describe children with relative coordinates, but allow computation of absolute ones
class ArticulatedBody {
	DynamicBody m_body;
	std::vector<pair<ArticulatedBody, Joint>> m_children;
};

#include "integration.hh"

// =====================

void on_key(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
		glfwSetWindowShouldClose(window, GL_TRUE);
		return;
	}
}

void on_mouse_button(GLFWwindow* window, int button, int action, int mods) {
	if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
	}
}

void on_scroll(GLFWwindow* window, double x, double y) {
}

// =====================
// OpenGL

int render_width = 0, render_height = 0;

void model_init(GLFWwindow* window) {
    glfwSetKeyCallback(window, on_key);
    glfwSetMouseButtonCallback(window, on_mouse_button);
    glfwSetScrollCallback(window, on_scroll);
}

Text* text = nullptr;

void render_init() {
	fprintf(stderr, "OpenGL version: [%s]\n", glGetString(GL_VERSION));
	glEnable(GL_CULL_FACE);
	glClearColor(0.2, 0.2, 1, 1.0);
	glViewport(0, 0, render_width, render_height);

	text = new Text;
	text->fg_color = vec4(1, 1, 1, 1);
	text->bg_color = vec4(0, 0, 0, 1);
}

void render_gui() {
	glm::mat4 matrix = glm::ortho<float>(0, render_width, 0, render_height, -1, 1);
	text->Reset(render_width, render_height, matrix, true);
	text->Print("Hello world!");
}

void OnError(int error, const char* message) {
	fprintf(stderr, "GLFW error %d: %s\n", error, message);
}

GLFWwindow* create_window() {
	glfwSetErrorCallback(OnError);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
	GLFWwindow* window = glfwCreateWindow(mode->width*2, mode->height*2, "Sima", glfwGetPrimaryMonitor(), NULL);
	if (!window)
		return nullptr;
	glfwMakeContextCurrent(window);
	glfwSwapInterval(0/*VSYNC*/);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwGetFramebufferSize(window, &render_width, &render_height);
	return window;
}

void sigsegv_handler(int sig) {
 	fprintf(stderr, "Error: signal %d:\n", sig);
    void* array[20];
    backtrace_symbols_fd(array, backtrace(array, 20), STDERR_FILENO);
    exit(1);
}

int main(int argc, char** argv) {
	void sigsegv_handler(int sig);
	signal(SIGSEGV, sigsegv_handler);
	if (!glfwInit())
        return -1;

	glm::dvec3 a;
	glm::i64vec3 b;
	a.x = glfwGetTime();
	b.x = rdtsc();
	usleep(100000);
	a.y = glfwGetTime();
	b.y = rdtsc();
	usleep(100000);
	a.z = glfwGetTime();
	b.z = rdtsc();

	GLFWwindow* window = create_window();
	if (!window)
        return -1;
	model_init(window);
	render_init();

	glm::dvec3 c;
	glm::i64vec3 d;
	c.x = glfwGetTime();
	d.x = rdtsc();
	usleep(100000);
	c.y = glfwGetTime();
	d.y = rdtsc();
	usleep(100000);
	c.z = glfwGetTime();
	d.z = rdtsc();
	Timestamp::init(a, b, c, d);

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		//model_frame(window, frame_ms);
		/*glm::mat4 matrix = glm::translate(perspective_rotation, -g_player.position);
		Frustum frustum(matrix);
		g_player.cpos = glm::ivec3(glm::floor(g_player.position)) >> ChunkSizeBits;*/

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);
		//render_world_blocks(matrix, frustum);
		glDisable(GL_DEPTH_TEST);
		render_gui();

		glfwSwapBuffers(window);
	}

	glfwTerminate();
	_exit(0); // exit(0) is not enough
	return 0;
}
