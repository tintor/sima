#include "generators.h"
#include "tesselate.h"
#include "align_alloc.h"
#include "auto.h"

mesh3 generate_box(double4 size) {
	aligned_vector<double4> vertices;
	vertices.reserve(8);
	for (double x : {-1, 1})
		for (double y : {-1, 1})
			for (double z : {-1, 1})
				vertices.push_back(double4{x, y, z, 1} * size);
	return convex_hull(vertices);
}

mesh3 generate_cylinder(uint sides, double rmin, double rmax, double zmin, double zmax) {
	aligned_vector<double4> vertices;
	for (auto [z, r] : {pair{zmin, rmin}, pair{zmax, rmax}})
		if (r == 0)
			vertices.push_back(double4{0, 0, z, 1});
		else
			for (int i : range(sides)) {
				double a = (2 * M_PI / sides) * i;
				double x = cos(a) * r;
				double y = sin(a) * r;
				vertices.push_back(double4{x, y, z, 1});
			}
	return convex_hull(vertices);
}

void AddQuad(mesh3& mesh, double4 a, double4 b, double4 c, double4 d) {
	mesh.emplace_back(a, b, c);
	mesh.emplace_back(c, d, a);
}

mesh3 CreateCrossMesh(double inner, double outer) {
	mesh3 m;
	double a = inner, b = outer;
	// TODO finish -> start with cube and use extrude on every face
	return m;
}

mesh3 generate_prism(const polygon2& poly, double zmin, double zmax) {
	mesh2 m2 = tesselate(poly);
	mesh3 m3;
	m3.reserve(m2.size() * 2 + poly.size() * 2);
	for (triangle2 m : m2) {
		// TODO check orientation of m
		m3.emplace_back(double4{m.a.x, m.a.y, zmin, 1}, double4{m.b.x, m.b.y, zmin, 1}, double4{m.c.x, m.c.y, zmin, 1});
		m3.emplace_back(double4{m.b.x, m.b.y, zmax, 1}, double4{m.a.x, m.a.y, zmax, 1}, double4{m.c.x, m.c.y, zmax, 1});
	}
	double2 sa = poly.back();
	for (double2 sb : poly) {
		double4 a{sa.x, sa.y, zmin, 1};
		double4 b{sa.x, sa.y, zmax, 1};
		double4 c{sb.x, sb.y, zmax, 1};
		double4 d{sb.x, sb.y, zmin, 1};
		m3.emplace_back(a, b, c);
		m3.emplace_back(c, d, a);
		sa = sb;
	}
	return m3;
}

double triangle_volume(double4 a, double4 b, double4 c);

struct Edge {
	int a, b;
	double len;
};

void print_volume(int e, span<const double4> vertex, span<const Edge> edges) {
	double volume = 0;
	for (auto i : range(edges.size() / 3)) {
		auto a = vertex[edges[i * 3].a];
		auto b = vertex[edges[i * 3 + 1].a];
		auto c = vertex[edges[i * 3 + 2].a];
		volume += triangle_volume(a, b, c);
	}
	volume /= 6;
	print("%s: volume %s\n", e, abs(volume));
}

mesh3 generate_regular_polyhedra(span<const Edge> edges) {
	int count = 0;
	for (auto e : edges)
		count = max(count, e.a + 1, e.b + 1);

	AlignAlloc<double4, 32> mem;
	double4* vertex = mem.allocate(count);
	ON_SCOPE_EXIT(mem.deallocate(vertex, count));
	double4* delta = mem.allocate(count);
	ON_SCOPE_EXIT(mem.deallocate(delta, count));

	std::default_random_engine rnd;
	rnd.seed(time(nullptr));
	for (auto i : range(count))
		vertex[i] = uniform3(rnd, -1, 1);

	print_volume(0, span<const double4>(vertex, count), edges);
	for (auto e : range(100)) {
		for (auto i : range(count))
			delta[i] = {0, 0, 0};
		for (auto e : edges) {
			double4 d = vertex[e.a] - vertex[e.b];
			d *= 0.05 * (e.len / length(d) - 1);
			delta[e.a] += d;
			delta[e.b] -= d;
		}
		for (auto i : range(count))
			vertex[i] += delta[i];
		print_volume(e + 1, span<const double4>(vertex, count), edges);
	}

	return mesh3();
}

inline double2 vec(int k, int n) {
	double a = (2 * M_PI * k) / n;
	return {cos(a), sin(a)};
}

mesh3 generate_regular_polyhedra(span<const span<const int>> faces) {
	vector<Edge> edges;
	for (auto f : faces) {
		int c = 1;
		double2 cv = vec(c, f.size());
		for (int a = 2; a < f.size(); a++) {
			int b = (a + 1) % f.size();
			double2 av = vec(a, f.size());
			double2 bv = vec(b, f.size());
			double len = length(av - bv);
			edges.push_back({f[a], f[b], 1.0});
			edges.push_back({f[b], f[c], length(bv - cv) / len});
			edges.push_back({f[c], f[a], length(cv - av) / len});
		}
	}
	return generate_regular_polyhedra(edges);
}

double print_volume2(int e, span<const double4> vertex, unordered_map<pair<int, int>, double>& edges) {
	double loss = 0;
	const double Inf = std::numeric_limits<double>::max();
	// for each vertex find its top3 nearest vertices (#1 and #2 should be at distance 0, #3 at distance 1)
	for (auto i : range(vertex.size())) {
		double dist[3] = {Inf, Inf, Inf};
		for (auto j : range(vertex.size()))
			if (j != i) {
				double d = length(vertex[i] - vertex[j]);
				if (d < dist[2]) {
					dist[2] = d;
					if (d < dist[1]) {
						swap(dist[1], dist[2]);
						if (d < dist[0])
							swap(dist[0], dist[1]);
					}
				}
			}
		if (e % 10000 == 9999)
			print("vertex[%s] = (%s), %s %s %s\n", i, vertex[i], dist[0], dist[1], dist[2]);
		loss += squared(dist[0]) + squared(dist[1]) + squared(1 - dist[2]);
	}
	if (e % 100 == 99)
		print("iter %s - %s\n", e, loss);
	return loss;
}

// create each regular polygon face first as free floating and all diagonal edges
// create weaker attraction force between every pair of edges (not on the same face)
// if not sealed after 100 iterations, randomize and repeat
mesh3 generate_regular_polyhedra2(span<const pair<int, int>> faces) {
	int count = 0;
	for (auto p : faces)
		count += p.first * p.second;

	aligned_vector<double4> vertex(count), delta(count);

	unordered_map<pair<int, int>, double> edges;
	count = 0;
	for (auto f : faces) {
		double len = length(vec(0, f.first) - vec(1, f.first));
		for (auto ii : range(f.second)) {
			// add a new face
			for (auto a : range(f.first))
				for (auto b : range(a)) {
					double2 av = vec(a, f.first), bv = vec(b, f.first);
					edges[pair<int, int>(count + a, count + b)] = length(av - bv) / len;
				}
			count += f.first;
		}
	}

	std::default_random_engine rnd;
	rnd.seed(time(nullptr));
	while(true) {
		// TODO simulate individual faces are rigid bodies, with attractive force between edges of different faces
		// TODO better init face vertices (faces away from each other)
		for (auto i : range(count))
			vertex[i] = uniform3(rnd, -1, 1);

		print_volume2(0, vertex, edges);
		for (auto e : range(200000)) {
			for (auto i : range(count))
				delta[i] = {0, 0, 0};
			for (auto a : range(count))
				for (auto b : range(a)) {
					auto it = edges.find(pair(a, b));
					double4 d = vertex[a] - vertex[b];
					double dlen = length(d);
					d /= dlen;
					if (it != edges.end()) {
						d *= 0.05 * (it->second - dlen);
					} else {
						d *= -min(0.5, 0.00001 / (dlen * dlen));
					}
					delta[a] += d;
					delta[b] -= d;
				}
			for (auto i : range(count))
				vertex[i] += delta[i];
			print_volume2(e + 1, vertex, edges);
		}
		break;
	}
	return mesh3();
}
