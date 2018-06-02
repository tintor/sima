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

//#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "tinyformat.h"
#include "auto.h"

using std::pair;
using std::cerr;
using std::cout;
using std::endl;
using std::size_t;
using std::min;
using std::max;

#include "primitives.hh"

typedef std::vector<triangle3> Mesh3d;

// Valid if triangles are not intersecting, except in one shared edge or one shared vertex
bool AreValidMeshFaces(const triangle3& a, const triangle3& b) {
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
			if (squared(*a - *b) <= 100 * squared(PlanarEpsilon))
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
			if (!AreValidMeshFaces(mesh[i], mesh[j]))
				return Validity::SelfIntersection;

	return Validity::OK;
}

TEST_CASE("is_valid") {
	REQUIRE(is_valid(std::vector<triangle3>{}) == Validity::TooFewFaces);
	dvec3 a(0, 0, 0);
	dvec3 b(1, 0, 0);
	dvec3 c(0, 1, 0);
	dvec3 d(0, 0, 1);
}

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
		// TODO can use the same trick from ConvexHull to avoid computing 1/sqrt()
		plane plane(f);
		FOR_EACH(v, vertices)
			if (plane.distance(v) > PlanarEpsilon)
				return false;
	}
	return true;
}

// Mesh generators
// ===============

// Returns empty vector if no solution (points are coplanar)
Mesh3d ConvexHullMesh(const std::vector<dvec3>& points) {
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
	return ConvexHullMesh(vertices);
}

Mesh3d SphereMesh(int vertices) {
	std::default_random_engine random;
	std::normal_distribution<real> gauss(0.0, 1.0);
	std::vector<dvec3> V(vertices);
	FOR(i, vertices)
		V[i] = glm::normalize(dvec3(gauss(random), gauss(random), gauss(random)));
	return ConvexHullMesh(V);
}

TEST_CASE("SphereMesh - ConvexHull - IsConvex") {
	Mesh3d mesh = SphereMesh(100);
	real v = 4.0 / 3 * M_PI;
	REQUIRE(abs(volume(mesh) - v) / v < 0.15);
	REQUIRE(is_valid(mesh) == Validity::OK);
	REQUIRE(is_convex(mesh));
}

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

TEST_CASE("CreateCrossMesh") {
	Mesh3d cross = CreateCrossMesh(0.05, 0.25);
	//REQUIRE(IsValid(cross) == Validity::OK);
}

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
		m_box_min = m_box_max = dvec3(0, 0, 0);
		FOR_EACH(v, m_convex_vertices) {
			m_sphere_radius = std::max(m_sphere_radius, squared(v));
			m_box_min = min(m_box_min, v);
			m_box_max = max(m_box_max, v);
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
	real sphere_radius() const { return m_sphere_radius; }

private:
	// TODO move faces / edges / vertices that collide more often to the start of the vector
	Mesh3d m_mesh; // TODO to save memory: keep unique vertices in separate array and use 16bit pointers
	std::vector<segment3> m_convex_edges; // non-concave and non-planar
	std::vector<dvec3> m_convex_vertices; // non-concave and non-planar (needs to have more than a hemisphere around it open)

	dmat3 m_inertia_tensor; // assuming density of 1kg/m^3
	real m_volume;
	real m_sphere_radius;
	dvec3 m_box_min, m_box_max; // object oriented bounding box

	bool m_is_convex;
	// TODO std::vector<plane> m_convex_planes;
	// TODO optimize if simple box mesh
	// TODO maybe plane for each face? (or at least a unit normal)
	// TODO maybe bounding box: AABB or OBB?
	// TODO maybe reduced bounding convex hull that fully encloses mesh?
	// TODO BSP tree for faster vertex VS mesh test
	// TODO voxel tensor for fastedr vertex VS mesh test
};

const auto random_directions = [](){
	std::array<dvec3, 128> dirs;
	std::default_random_engine rnd;
	FOR_EACH(d, dirs)
		d = random_direction(rnd);
	return dirs;
}();

struct Body {
	Shape shape;
	dvec3 position;
	dquat orientation;
};

bool IsVertexPenetratingShape(const dvec3& v, const Shape& shape) {
	// TODO bounding box / sphere check

	// TODO optimize for box
	// TODO optimize for convex mesh

	// TODO even if mesh is concave we can have a convex bounding shape to quickly test against
	// (only if convex hull has small number of faces)

	// TODO if possible to split concave mesh into small number of convex ones: do it
	//      (even splitting large concave into smaller concave helps)

	// TODO could be precomputed with a dense voxel table:
	// each voxel is one of: outside / convex maybe / concave maybe / inside
	// convex maybe has list of faces that intersect the voxel (1 face in best case)
	// in case of concave maybe run full ray intersect algorithm

	// TODO optimize convex planar faces with more than 3 vertices (ie. sides of a cube)

	FOR_EACH(dir, random_directions) {
		int hits = 0;
		FOR_EACH(face, shape.faces()) {
			// TODO TODO compare ray (v,dir) and face
			// if ray strictly goes through the triangle then hits ^= 1
			// else if ray intersects face plane in point outside of face then continue
			// else if ray is coplanar with face but far from face bounding circle then continue
			// else (ray hits edge, vertex or is coplanar with face) then restart while loop
		}
		return hits != 0;
	}
    throw new std::runtime_error("IsVertexPenetratingMesh failure");
}

bool IsEdgePenetratingShape(const segment3& edge, const Shape& shape) {
	// TODO TODO At least one of these negative cases is needed for performance
	// TODO use voxel grid here faster positive and negative cases

	// TODO bounding box / sphere check
	// TODO optimize for box
	// TODO optimize for convex mesh

	real edge_length = l2Norm(edge.b - edge.a);
	// t_min and t_max needed to avoid edge endpoints being too close triangle an merely touching it
	real t_min = ContactEpsilon / edge_length;
	real t_max = 1 - t_min;

	FOR_EACH(face, shape.faces()) {
		// compare edge and triangle
		dvec3 normal = Normal(face);
		real dd = dot(normal, edge.b - edge.a);
		// if edge is parallel to triangle plane
		constexpr real tiny = 1e-8;
		if (abs(dd) < tiny || abs(dd) < tiny * l2Norm(normal) * edge_length)
			continue;
		real t = dot(normal, face.a - edge.a) / dd;
		// if intersection of infinite edge with triangle plane is outside of edge interior
		if (t < t_min || t > t_max)
			continue;
		dvec3 w = edge.linear(t); // intersection of triangle plane and edge

		FOR_EACH_EDGE(a, b, face) {
			// if W is outside of triangle
			if (plane::sign(*a, *b, *a + normal, w) > 0)
				return false;
			// if edge is too close to triangle edge (could be contact instead)
			if (segment3::squared_distance(edge, segment3(*a, *b)) <= squared(ContactEpsilon))
				return false;
		}
		return true;
	}

	// Hard cases:
	// - edges don't go through interior of any triangles
	// - endpoints aren't inside polyhedron
	//   => endpoints can be outside and on the surface
	// segment can intersect surface at the vertex or interior of edge

	// Ignoring hard cases as they can't happen in simulator (only adversarial examples)

	// Idea: look at immediate surface around polyhedron vertex or edge that intersects segment
	// That surface can be convex
	return false;
}

struct OBB {
    dvec3 size;
};

pair<real, real> project_obb(dvec3 obb_position, dvec3 obb_size, dmat3 obb_orientation, dvec3 dir) {
    // TODO
	return pair<real, real>(0.0, 0.0);
}

bool intersects_obb(const Body& body_a, const Body& body_b) {
    // TODO
	return false;
}

bool AreShapesSeparatedQuick(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
    // TODO sphere / sphere
    // TODO OBB / OBB
    // TODO convex / convex (if at least one of them is concave)
	return false;
}

bool AreBodiesPenetrating(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	// Check all edges of one against the other
	FOR_EACH(edge_a, shape_a.convex_edges())
		if (IsEdgePenetratingShape(pos_b.to_local(pos_a.to_global(edge_a)), shape_b))
				return true;
	FOR_EACH(edge_b, shape_b.convex_edges())
		if (IsEdgePenetratingShape(pos_a.to_local(pos_b.to_global(edge_b)), shape_a))
				return true;

	// Also check a single vertex from one body against another,
	// in case of one shape being completely inside the other.
	return IsVertexPenetratingShape(pos_b.to_local(pos_a.to_global(shape_a.convex_vertices()[0])), shape_b)
		|| IsVertexPenetratingShape(pos_a.to_local(pos_b.to_global(shape_b.convex_vertices()[0])), shape_a);
}

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

// Assuming bodies aren't penetrating each other
std::vector<Contact> find_all_contacts(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	std::vector<Contact> contacts;
	Contact contact;

	// vertex vs face contacts
	FOR_EACH(vertex_a, shape_a.convex_vertices()) {
		dvec3 vertex_a_global = pos_a.to_global(vertex_a);
		// TODO bounding volume check?
		dvec3 vertex_a_local = pos_b.to_local(vertex_a_global);
		FOR_EACH(face_b, shape_b.faces())
			if (is_vertex_triangle_contact(vertex_a_local, face_b, /*out*/contact)) {
				contact.position = pos_b.to_global(contact.position);
				contact.normal = pos_b.to_global_dir(contact.normal);
				contacts.push_back(contact);
			}
	}
	FOR_EACH(vertex_b, shape_b.convex_vertices()) {
		dvec3 vertex_b_global = pos_b.to_global(vertex_b);
		// TODO bounding volume check?
		dvec3 vertex_b_local = pos_a.to_local(vertex_b_global);
		FOR_EACH(face_a, shape_a.faces())
			if (is_vertex_triangle_contact(vertex_b_local, face_a, /*out*/contact)) {
				contact.position = pos_a.to_global(contact.position);
				contact.normal = -pos_a.to_global_dir(contact.normal);
				contacts.push_back(contact);
			}
	}

	// edge vs edge contacts
	FOR_EACH(edge_a, shape_a.convex_edges()) {
		segment3 edge_a_global = pos_a.to_global(edge_a);
		// TODO bounding volume check?
		// TODO avoid transforming every edge in B
		FOR_EACH(edge_b, shape_b.convex_edges())
			if (is_edge_edge_contact(edge_a_global, pos_b.to_global(edge_b), /*out*/contact))
				contacts.push_back(contact);
	}

	return contacts;
}

bool are_bounding_spheres_intersecting(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
	return squared(pos_a.position - pos_b.position) <= squared(shape_a.sphere_radius() + shape_b.sphere_radius());
}

// returns 0 if interecting
real distance(const Shape& shape_a, const transform3& pos_a, const Shape& shape_b, const transform3& pos_b) {
    // TODO need strict intersection (not just deep penetration)
    if (!are_bounding_spheres_intersecting(shape_a, pos_a, shape_b, pos_b)
			&& IsVertexPenetratingShape(pos_b.to_local(pos_a.to_global(shape_a.convex_vertices()[0])), shape_b))
        return 0;

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
