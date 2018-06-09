#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <regex>

#include "shape.hh"
#include "util.hh"
#include "auto.h"
#include "file.hh"

constexpr const char* FLOAT_REGEX = R"(-?\d+\.\d+(e[+-]\d+)?)";

Mesh3d load_stl(const std::string& filename) {
    FileReader file(filename.c_str());
    std::string_view line;

    Mesh3d mesh;
    std::cmatch match;

    if (!file.getline(line) || line != "solid Default")
        throw std::runtime_error("expected [solid Default]");

    std::regex facet_regex(R"(^\s*facet\s+)");
    std::regex vertex_regex(tfm::format(R"(^\s*vertex (%s) (%s) (%s)\s*$)", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (true) {
        if (!file.getline(line))
            throw std::runtime_error("expected [facet ...] or [endsolid Default] insted of end");
        if (line == "endsolid Default")
            break;
        if (!std::regex_search(line.begin(), line.end(), facet_regex))
            throw std::runtime_error(tfm::format("expected [facet ...] or [endsolid Default] instead of [%s]", line));
        if (!file.getline(line) || line != "    outer loop")
            throw std::runtime_error("expected [outer loop]");

        dvec3 v[3];
        FOR(i, 3) {
            if (!file.getline(line) || !std::regex_match(line.begin(), line.end(), match, vertex_regex))
                throw std::runtime_error(tfm::format("expected [vertex ...] instead of [%s]", line));

            from_chars(match[1].first, match[1].second, v[i].x);
            from_chars(match[3].first, match[3].second, v[i].y);
            from_chars(match[5].first, match[5].second, v[i].z);
        }
        mesh.emplace_back(v[0], v[1], v[2]);

        if (!file.getline(line) || line != "    endloop")
            throw std::runtime_error("expected [endloop]");
        if (!file.getline(line) || line != "  endfacet")
            throw std::runtime_error("expected [endfacet]");
    }

    return mesh;
}

Mesh3d load_ply(const std::string& filename) {
    FileReader file(filename.c_str());
    std::string_view line;

    std::vector<dvec3> vertices;
    Mesh3d mesh;
    std::cmatch match;

    std::regex element_regex(R"(element (vertex|face) (\d+)\s*)");
    while (true) {
        if (!file.getline(line))
            throw std::runtime_error("bad ply file header");
        if (line == "end_header")
            break;
        if (std::regex_match(line.begin(), line.end(), /*out*/match, element_regex)) {
            int s;
            from_chars(match[2].first, match[2].second, s);
            if (match[1].str() == "vertex")
                vertices.reserve(s);
            else
                mesh.reserve(s);
        }
    }

    std::regex vertex_regex(tfm::format(R"(^(%s) (%s) (%s)(\s+|$))", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (vertices.size() < vertices.capacity()) {
        if (!file.getline(line))
            throw std::runtime_error("unexpected end of ply file");
        if (!std::regex_search(line.begin(), line.end(), /*out*/match, vertex_regex))
            throw std::runtime_error(tfm::format("expected vertex in ply file: [%s]", line));
        double a, b, c;
        from_chars(match[1].first, match[1].second, a);
        from_chars(match[3].first, match[3].second, b);
        from_chars(match[5].first, match[5].second, c);
        vertices.emplace_back(a, b, c);
    }

    std::regex face_regex(R"(3 (\d+) (\d+) (\d+)\s*)");
    while (mesh.size() < mesh.capacity()) {
        if (!file.getline(line) || !std::regex_match(line.begin(), line.end(), /*out*/match, face_regex))
            throw std::runtime_error("expected face in ply file");
        int a, b, c;
        from_chars(match[1].first, match[1].second, a);
        from_chars(match[3].first, match[3].second, b);
        from_chars(match[5].first, match[5].second, c);
        mesh.emplace_back(vertices[a], vertices[b], vertices[c]);
    }

    return mesh;
}

// --------------------

const SolidLeafBSPTree* SolidLeafBSPTree::Inside = reinterpret_cast<const SolidLeafBSPTree*>(1);

int classify(dvec3 v, const plane& p, double epsilon) {
	auto d = p.distance(v);
	if (d > epsilon)
		return 1;
	if (d < -epsilon)
		return -1;
	return 0;
}

enum class TriangleClass : uint8_t {
	Positive,
	Negative,
	Coplanar,
	Split,
};

TriangleClass classify(const triangle3& f, const plane& p, double epsilon) {
    int a = classify(f.a, p, epsilon);
    int b = classify(f.b, p, epsilon);
    int c = classify(f.c, p, epsilon);

    int mn = min(a, b, c);
    int mx = max(a, b, c);

    if (mx == 1 && mn != -1)
        return TriangleClass::Positive;
    if (mn == -1 && mx != 1)
        return TriangleClass::Negative;
    if (mn == 0 && mx == 0)
        return TriangleClass::Coplanar;
    return TriangleClass::Split;
}

struct SolidBSPTreeBuildData {
	std::vector<dvec3> vertices;
	std::vector<dvec3> samples;
	std::default_random_engine rnd;
	std::vector<TriangleClass> side;
};

int depth = 0;

// TODO add memoization / add persisted memoization
// Returns pair of tree node and average query depth across all samples
std::pair<SolidLeafBSPTree*, real> build_solid_leaf_bsp_tree_internal(
        SolidBSPTreeBuildData& data, const std::vector<uint16_t>& mesh,
        std::vector<uint32_t>::iterator samples_begin,
        std::vector<uint32_t>::iterator samples_end) {
    depth += 1;
    AUTO(depth -= 1);
    std::cout << "build " << (mesh.size() / 3) << " (" << depth << ")" << std::endl;
	if (mesh.empty())
		return std::pair<SolidLeafBSPTree*, real>(nullptr, 0.0);

    if (mesh.size() == 3) {
        SolidLeafBSPTree* m = new SolidLeafBSPTree;
        m->divider = plane(data.vertices[mesh[0]], data.vertices[mesh[1]], data.vertices[mesh[2]]);
        m->positive = nullptr;  // outside
        m->negative = const_cast<SolidLeafBSPTree*>(SolidLeafBSPTree::Inside);
        return std::pair<SolidLeafBSPTree*, real>(m, 1.0);
    }

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
    release(vertices);

	SolidLeafBSPTree* tree = new SolidLeafBSPTree;
	real best_score = 1e100;
	std::vector<uint16_t> pmesh, nmesh;

    std::uniform_int_distribution<uint32_t> uint_dist(0, candidates.size() - 1);
    //FOR_EACH(candidate, candidates)
    int c = uint_dist(data.rnd);
    while (true) {
        auto candidate = candidates[c++ % candidates.size()];
		// classify faces against divider and build tree recursively
        data.side.resize(mesh.size() / 3);
        int s = mesh.size() / 3;
        FOR(f, mesh.size() / 3) {
            data.side[f] = classify(triangle3(data.vertices[mesh[f * 3]], data.vertices[mesh[f * 3 + 1]], data.vertices[mesh[f * 3 + 2]]), candidate, 1e-5);
            if (data.side[f] == TriangleClass::Coplanar)
                s -= 1;
            if (data.side[f] == TriangleClass::Split)
                s += 1;
        }

        if (s <= mesh.size() / 3) {
            // repartition in place
        } else {
            // use pmesh and nmesh
        }

        pmesh.clear();
		nmesh.clear();
        FOR(f, mesh.size() / 3) {
			if (data.side[f] == TriangleClass::Positive || data.side[f] == TriangleClass::Split)
				FOR(i, 3)
					pmesh.push_back(mesh[f * 3 + i]);
            if (data.side[f] == TriangleClass::Negative || data.side[f] == TriangleClass::Split)
				FOR(i, 3)
					nmesh.push_back(mesh[f * 3 + i]);
		}

        if (pmesh.size() == 0 || nmesh.size() == 0)
            if (pmesh.size() + nmesh.size() == mesh.size())
                continue;
        std::cout << (pmesh.size() / 3) << " : " << (nmesh.size() / 3) << std::endl;

        auto split = std::partition(samples_begin, samples_end, [&data, &candidate](uint32_t s) {
            return candidate.distance(data.samples[s]) > 0;
        });
		auto presult = build_solid_leaf_bsp_tree_internal(data, pmesh, samples_begin, split);
		auto nresult = build_solid_leaf_bsp_tree_internal(data, nmesh, split, samples_end);

        auto p = std::distance(samples_begin, split) * presult.second;
        auto n = std::distance(split, samples_end) * nresult.second;
		real score = 1 + (p + n) / std::distance(samples_begin, samples_end);

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
        break;
	}
	return std::pair(tree, best_score);
}

dvec3 min(const Mesh3d& mesh) {
    dvec3 vmin = mesh[0].a;
    FOR_EACH(face, mesh)
		FOR(i, 3)
            vmin = min(vmin, face[i]);
	return vmin;
}

dvec3 max(const Mesh3d& mesh) {
    dvec3 vmax = mesh[0].a;
    FOR_EACH(face, mesh)
		FOR(i, 3)
            vmax = max(vmax, face[i]);
	return vmax;
}

std::unique_ptr<SolidLeafBSPTree> build_solid_leaf_bsp_tree(const Mesh3d& mesh, uint32_t num_samples) {
	SolidBSPTreeBuildData data;

    dvec3 vmin = min(mesh);
    dvec3 vmax = max(mesh);
    dvec3 d = vmax - vmin;
    real dmax = max(d.x, d.y, d.z);

	// Init samples
	data.samples.resize(num_samples);
	std::vector<uint32_t> isamples(num_samples);
	FOR(i, num_samples) {
        while (true) {
            dvec3 v = random_vector(data.rnd) * d + vmin;
            if (v.x <= vmax.x && v.y <= vmax.y && v.z <= vmax.z) {
                data.samples[i] = v;
        		isamples[i] = i;
                break;
            }
        }
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

	auto result = build_solid_leaf_bsp_tree_internal(data, imesh, isamples.begin(), isamples.end());
    std::cout << "TreeScore " << result.second << std::endl;
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

bool is_box(const Mesh3d& mesh) {
	dvec3 vmin = min(mesh);
	dvec3 vmax = max(mesh);

	// All vertices must be made from extreme coordinates
	FOR_EACH(f, mesh)
		FOR(i, 3) {
			dvec3 v = f[i];
			if (v.x != vmin.x && v.x != vmax.x && v.y != vmin.y && v.y != vmax.y && v.z != vmin.z && v.z != vmax.z)
				return false;
		}

	// Every face must have one coordinate constant
	FOR_EACH(f, mesh) {
		bool xx = f.a.x == f.b.x && f.b.x == f.c.x;
		bool yy = f.a.y == f.b.y && f.b.y == f.c.y;
		bool zz = f.a.z == f.b.z && f.b.z == f.c.z;
		if (!xx && !yy && !zz)
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

// Shape is 3d solid, immutable, purely geometric and with origin in center of mass
Shape::Shape(const Mesh3d& mesh, /*in/out*/transform3& position) {
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

	// Bounding box
	m_is_box = ::is_box(m_mesh);
	m_box.first = min(m_mesh);
	m_box.second = max(m_mesh);

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

	// Compute radius of bounding sphere
	m_sphere_radius = 0;
	FOR_EACH(v, m_convex_vertices)
		m_sphere_radius = std::max(m_sphere_radius, squared(v));
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

/*const auto random_directions = [](){
	std::array<dvec3, 128> dirs;
	std::default_random_engine rnd;
	FOR_EACH(d, dirs)
		d = random_direction(rnd);
	return dirs;
}();*/

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

std::pair<real, real> project_obb(dvec3 obb_position, dvec3 obb_size, dmat3 obb_orientation, dvec3 dir) {
    // TODO
	return std::pair<real, real>(0.0, 0.0);
}

bool are_oriented_boxes_intersecting(const std::pair<dvec3, dvec3>& box_a, const transform3& pos_a, const std::pair<dvec3, dvec3>& box_b, const transform3& pos_b) {
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
