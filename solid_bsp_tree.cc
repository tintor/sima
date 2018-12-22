#include "solid_bsp_tree.h"
#include "aabb.h"
#include <core/util.h>
#include "properties.h"
#include <core/exception.h>
#include "primitives.h"

void SolidBSPTree::evaluate_candidate(
		plane candidate,
		BuildData& data,
		const vector<uint>& mesh,
		uint* samples_begin,
		uint* samples_end,
		long& best_heuristic,
		Hist& best_hist,
		plane& best_candidate) {
	size_t psamples = 0;
	for (auto it = samples_begin; it != samples_end; it++)
		if (candidate.distance(data.samples[*it]) > 0)
			psamples += 1;
	size_t samples = samples_end - samples_begin;
	size_t nsamples = samples - psamples;
	if (psamples == 0 || nsamples == 0)
		return;

	Hist hist{ 0, 0, 0, 0 };
	long h = best_heuristic;
	for (auto f : mesh) {
		const triangle3& face = data.faces[mesh[f]];

		int ai = candidate.classify(face.a);
		int bi = candidate.classify(face.b);
		int ci = candidate.classify(face.c);

		int mn = min(ai, bi, ci);
		int mx = max(ai, bi, ci);

		if (mx == 1 && mn != -1) {
			hist.positive += 1;
			h -= psamples;
			if (h <= 0)
				return;
			continue;
		}
		if (mn == -1 && mx != 1) {
			hist.negative += 1;
			h -= nsamples;
			if (h <= 0)
				return;
			continue;
		}
		if (mn == -1 && mx == 1) {
			hist.stradle += 1;
			h -= samples;
			if (h <= 0)
				return;
			continue;
		}
	}

	if (h > 0) {
		// heuristic = psamples * hist.positive + nsamples * hist.negative + samples * hist.stradle;
		best_heuristic -= h;
		hist.overlap = mesh.size() - hist.positive - hist.negative - hist.stradle;
		best_hist = hist;
		best_candidate = candidate;
	}
}

// Returns average query depth across all samples
pair<uint, float> SolidBSPTree::build_internal(
		BuildData& data, const vector<uint>& mesh, uint* samples_begin, uint* samples_end) {
	const uint Outside = 0;
	const uint Inside = 0xFFFFFFFF;

	if (mesh.size() >= 1000)
		print("build_internal: %s\n", mesh.size());

	if (mesh.empty())
		return pair<uint, float>(Outside, 0.0);

	if (mesh.size() == 1) {
		int i = add_node();
		node[i].divider = plane(data.faces[mesh[0]]);
		node[i].positive = Outside;
		node[i].negative = Inside;
		percent[i] = float(std::distance(samples_begin, samples_end)) / data.samples.size();
		hist[i] = { 0, 0, 0, 0 };
		hist[i].overlap = 1;
		box_size[i] = aabb4(data.faces[mesh[0]]).size();
		return pair<uint, double>(i, 1.0);
	}

	vector<uint16_t> vertices;
	vertices.reserve(mesh.size() * 3);
	for (auto f : mesh) {
		vertices.push_back(data.ifaces[f].a);
		vertices.push_back(data.ifaces[f].b);
		vertices.push_back(data.ifaces[f].c);
	}
	remove_dups(vertices);

	const auto samples = std::distance(samples_begin, samples_end);
	long best_heuristic = std::numeric_limits<long>::max();
	Hist best_hist { 0, 0, 0, 0 };
	plane best_candidate;

	// TODO if Faces * Samples is too high then skip that too, and just use linear major axis algo

	// if V is too high then fallback to O(Vertices + Faces * Samples), use every face as candidate + 3 main planes through every vertex
	if (vertices.size() > 10) {
		// Consider only subset of candidates
		if (mesh.size() <= 10000) for (auto f : mesh) {
			// TODO precompute
			auto& v = data.vertices;
			plane candidate(vertices[data.ifaces[f].a], vertices[data.ifaces[f].b], vertices[data.ifaces[f].c]);
			evaluate_candidate(candidate, data, mesh, samples_begin, samples_end,
					best_heuristic, best_hist, best_candidate);
		}
		for (auto major_axis : range(3)) {
			double4 normal{0, 0, 0};
			normal[major_axis] = 1;

			// sort vertices along the major axis
			const double4* vertex_list = &(data.vertices[0]);
			std::sort(vertices.begin(), vertices.end(), [major_axis, vertex_list](uint16_t va, uint vb){
				return vertex_list[va][major_axis] < vertex_list[vb][major_axis];
			});

			// sort samples along the major axis
			const double4* sample_list = &(data.samples[0]);
			std::sort(samples_begin, samples_end, [major_axis, sample_list](uint sa, uint sb){
				return sample_list[sa][major_axis] < sample_list[sb][major_axis];
			});

			double v_prev = -1e100;
			for (auto v : vertices) {
				double4 vv = data.vertices[v];
				if (-vv[major_axis] == v_prev)
					continue;
				v_prev = vv[major_axis];
				plane candidate(normal, -vv[major_axis]);
				// TODO pre-sort samples outside of vertex loop to avoid re-partitioning
				evaluate_candidate(candidate, data, mesh, samples_begin, samples_end,
						best_heuristic, best_hist, best_candidate);
			}
		}
	} else {
		// Consider all candidates
		// Expensive! O(Vertices^3 * Samples)
		for (auto a : range(vertices.size())) for (auto b : range(a)) for (auto c : range(b)) {
			plane candidate(data.vertices[vertices[a]], data.vertices[vertices[b]], data.vertices[vertices[c]]);
			evaluate_candidate(candidate, data, mesh, samples_begin, samples_end,
					best_heuristic, best_hist, best_candidate);
		}
	}

	auto split = samples_begin;
	if (best_heuristic == std::numeric_limits<long>::max()) {
		if (mesh.size() > 250)
			print("out of samples: %s\n", mesh.size());
		// Out of samples -> split using new heuristic ignoring samples: abs(positive - negative) + 8 * splits
		for (auto a : range(vertices.size())) for (auto b : range(a)) for (auto c : range(b)) {
			plane candidate(data.vertices[vertices[a]], data.vertices[vertices[b]], data.vertices[vertices[c]]);

			Hist hist{ 0, 0, 0, 0 };
			for (auto f : range(mesh.size())) {
				const triangle3& face = data.faces[mesh[f]];

				int ai = candidate.classify(face.a);
				int bi = candidate.classify(face.b);
				int ci = candidate.classify(face.c);

				int mn = min(ai, bi, ci);
				int mx = max(ai, bi, ci);

				if (mx == 1 && mn != -1) {
					hist.positive += 1;
					continue;
				}
				if (mn == -1 && mx != 1) {
					hist.negative += 1;
					continue;
				}
				if (mn == -1 && mx == 1) {
					hist.stradle += 1;
					continue;
				}
			}

			long heuristic = abs(hist.positive - hist.negative) + 8 * hist.stradle;
			if (heuristic < best_heuristic) {
				best_heuristic = heuristic;
				hist.overlap = mesh.size() - hist.positive - hist.negative - hist.stradle;
				best_hist = hist;
				best_candidate = candidate;
			}
		}
		split = samples_end = samples_begin;
	} else {
		split = std::partition(samples_begin, samples_end, [&data, &best_candidate](uint s) {
			return best_candidate.distance(data.samples[s]) > 0;
		});
	}

	const auto psamples = std::distance(samples_begin, split);
	const auto nsamples = std::distance(split, samples_end);

	vector<uint> pmesh, nmesh;
	pmesh.reserve(best_hist.positive + best_hist.stradle);
	nmesh.reserve(best_hist.negative + best_hist.stradle);
	for (auto f : mesh) {
		const triangle3& face = data.faces[f];

		int ai = best_candidate.classify(face.a);
		int bi = best_candidate.classify(face.b);
		int ci = best_candidate.classify(face.c);

		if (ai > 0 || bi > 0 || ci > 0)
			pmesh.push_back(f);
		if (ai < 0 || bi < 0 || ci < 0)
			nmesh.push_back(f);
	}

	aabb4 box;
	for (auto iv : vertices)
		box.add(data.vertices[iv]);

	uint i = add_node();
	vertices.shrink_to_fit();

	auto presult = build_internal(data, pmesh, samples_begin, split);
	if (pmesh.size() == 0) {
		// TODO may need to correct the sign in some cases
		//result.first = (candidate vs face it was made from) ? SolidLeafBSPTree::Inside : nullptr;
	}
	auto nresult = build_internal(data, nmesh, split, samples_end);
	if (nmesh.size() == 0) {
		// TODO
	}

	percent[i] = static_cast<double>(samples) / data.samples.size();
	node[i].divider = best_candidate;
	hist[i] = best_hist;
	node[i].positive = presult.first;
	node[i].negative = nresult.first;
	double score = 1 + (psamples * presult.second + nsamples * nresult.second) / samples;
	box_size[i] = box.size();
	return pair<uint, double>(i, score);
}

void SolidBSPTree::print_tree(uint n, int depth) {
	for (auto i : range(depth))
		print("\t");
	const uint Inside = 0xFFFFFFFF;
	const uint Outside = 0;

	if (n == Outside && depth != 0) {
		print("outside\n");
		return;
	}
	if (n == Inside) {
		print("inside\n");
		return;
	}
	print("plane %.3f%% [%d %d %d %d] (%s)\n", 100 * percent[n],
		(int)hist[n].positive, (int)hist[n].negative, (int)hist[n].overlap, (int)hist[n].stradle, box_size[n]);
	print_tree(node[n].positive, depth + 1);
	print_tree(node[n].negative, depth + 1);
}

bool intersects_convex(double4 v, const mesh3& mesh) {
	for (auto face : mesh) {
		double4 normal = compute_normal(face);
		if (dot(normal, v - face.a) > 0)
			return false;
	}
	return true;
}

bool intersects(double4 v, const triangle3* mesh_begin, const triangle3* mesh_end) {
	// Can't just take the sdist sign of min_dist as two faces that share the edge
	// can have the same dist but different sgn(sdist)
	double min_dist = std::numeric_limits<double>::max(), max_sdist = 0;
	constexpr double eps = 1e-5; // ??? arbitrary, need to account for rounding errors
	// when computing distance(vertex, face)
	for (auto m = mesh_begin; m != mesh_end; m += 1) {
		const triangle3& face = *m;
		plane p(face); // TODO precompute this

		double dist, sdist;
		for (auto [a, b] : face.edges())
			if (plane::sign(a, b, a + p.normal(), v) > 0) {
				dist = distance(v, segment3{a, b});
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
	return max_sdist <= 0;
}

bool intersects(double4 v, const mesh3& mesh) {
	return ::intersects(v, &mesh[0], &mesh[0] + mesh.size());
}

bool SolidBSPTree::intersects(double4 v) {
	const uint Inside = 0xFFFFFFFF;
	const uint Outside = 0;

	uint i = 0;
	while (true) {
		const Node& n = node[i];
		i = (n.divider.distance(v) > 0) ? n.positive : n.negative;
		if (i == Outside)
			return false;
		if (i == Inside)
			return true;
	}
}

SolidBSPTree::SolidBSPTree(const mesh3& mesh, uint num_samples, std::default_random_engine& rnd) {
	BuildData data;
	data.rnd = rnd;
	aabb4 box(mesh);

	// Init samples
	data.samples.resize(num_samples);
	vector<uint> isamples(num_samples);
	for(auto i : range(num_samples)) {
		data.samples[i] = uniform3(data.rnd, box);
		isamples[i] = i;
	}

	// Init vertices and imesh
	unordered_map<double4, uint16_t, hash_t<double4>, equal_t<double4>> vec_index;
	vector<uint> imesh;
	for (const auto& face : mesh) {
		const auto& a = face[0];
		if (vec_index.count(a) == 0) {
			vec_index[a] = data.vertices.size();
			data.vertices.push_back(a);
		}
		const auto& b = face[1];
		if (vec_index.count(b) == 0) {
			vec_index[b] = data.vertices.size();
			data.vertices.push_back(b);
		}
		const auto& c = face[2];
		if (vec_index.count(c) == 0) {
			vec_index[c] = data.vertices.size();
			data.vertices.push_back(c);
		}
		itriangle t;
		t.a = vec_index[a];
		t.b = vec_index[b];
		t.c = vec_index[c];
		data.ifaces.push_back(t);
		data.faces.push_back(face);
		imesh.push_back(imesh.size());
	}

	auto result = build_internal(data, imesh, &isamples[0], &isamples[0] + isamples.size());
	print("TreeScore %s ", result.second);
	print_tree(0);

	int ab = 0, aa = 0, bb = 0;
	double n = num_samples;
	vector<double4> tester;
	tester.resize(n);
	for (auto j : range(n))
		tester[j] = uniform3(data.rnd, box);

	for (auto j : range(n)) {
		auto v = tester[j];
		bool insideA = ::intersects(v, mesh);
		bool insideB = intersects(v);
		if (insideA)
			aa += 1;
		if (insideB)
			bb += 1;
		if (insideA == insideB)
			ab += 1;
	}

	double4 d = box.size();
	double v = Volume(mesh) / (d.x * d.y * d.z);
	print("Same %.05f, A %.05f, Tree %.05f Volume %.05f", ab / n, aa / n, bb / n, v);
}
