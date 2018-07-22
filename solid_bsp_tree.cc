#include "solid_bsp_tree.h"

/*int classify(dvec3 v, const plane& p, double epsilon) {
	auto d = p.distance(v);
	if (d > epsilon)
		return 1;
	if (d < -epsilon)
		return -1;
	return 0;
}

void SolidBSPTree::evaluate_candidate(
        const plane& candidate,
        BuildData& data,
        const std::vector<uint32_t>& mesh,
        uint32_t* samples_begin,
        uint32_t* samples_end,
        int64_t& best_heuristic,
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
    int64_t h = best_heuristic;
    FOR(f, mesh.size()) {
        const triangle3& face = data.faces[mesh[f]];

        constexpr real epsilon = 1e-5;
        int ai = classify(face.a, candidate, epsilon);
        int bi = classify(face.b, candidate, epsilon);
        int ci = classify(face.c, candidate, epsilon);

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

dvec3 min(const std::vector<dvec3>& vertices) {
    dvec3 vmin = vertices[0];
    FOR_EACH(v, vertices)
        vmin = min(vmin, v);
	return vmin;
}

dvec3 max(const std::vector<dvec3>& vertices) {
    dvec3 vmax = vertices[0];
    FOR_EACH(v, vertices)
        vmax = max(vmax, v);
	return vmax;
}

// Returns average query depth across all samples
std::pair<uint32, real> SolidBSPTree::build_internal(
        BuildData& data, const std::vector<uint32_t>& mesh, uint32_t* samples_begin, uint32_t* samples_end) {
    const uint32_t Outside = 0;
    const uint32_t Inside = 0xFFFFFFFF;

    if (mesh.size() >= 1000)
        std::cout << "build_internal: " << mesh.size() << std::endl;

	if (mesh.empty())
		return std::pair<uint32_t, real>(Outside, 0.0);

    if (mesh.size() == 1) {
        int i = add_node();
        node[i].divider = plane(data.faces[mesh[0]]);
        node[i].positive = Outside;
        node[i].negative = Inside;
        percent[i] = float(std::distance(samples_begin, samples_end)) / data.samples.size();
        hist[i] = { 0, 0, 0, 0 };
        hist[i].overlap = 1;
        box_size[i] = max(data.faces[mesh[0]]) - min(data.faces[mesh[0]]);
        return std::pair<uint32_t, real>(i, 1.0);
    }

	std::vector<uint16_t> vertices;
    vertices.reserve(mesh.size() * 3);
	FOR_EACH(f, mesh) {
		vertices.push_back(data.ifaces[f].a);
        vertices.push_back(data.ifaces[f].b);
        vertices.push_back(data.ifaces[f].c);
    }
    remove_dups(vertices);

    const auto samples = std::distance(samples_begin, samples_end);
    int64_t best_heuristic = std::numeric_limits<int64_t>::max();
    Hist best_hist { 0, 0, 0, 0 };
    plane best_candidate;

    // TODO if Faces * Samples is too high then skip that too, and just use linear major axis algo

    // if V is too high then fallback to O(Vertices + Faces * Samples), use every face as candidate + 3 main planes through every vertex
    if (vertices.size() > 10) {
        // Consider only subset of candidates
        if (mesh.size() <= 10000) FOR_EACH(f, mesh) {
            // TODO precompute
            plane candidate(data.vertices[data.ifaces[f].a], data.vertices[data.ifaces[f].b], data.vertices[data.ifaces[f].c]);
            evaluate_candidate(candidate, data, mesh, samples_begin, samples_end,
                    best_heuristic, best_hist, best_candidate);
        }
        FOR(major_axis, 3) {
            dvec3 normal(0, 0, 0);
            normal[major_axis] = 1;

            // sort vertices along the major axis
            const dvec3* vertex_list = &(data.vertices[0]);
            std::sort(vertices.begin(), vertices.end(), [major_axis, vertex_list](uint16_t va, uint32_t vb){
                return vertex_list[va][major_axis] < vertex_list[vb][major_axis];
            });

            // sort samples along the major axis
            const dvec3* sample_list = &(data.samples[0]);
            std::sort(samples_begin, samples_end, [major_axis, sample_list](uint32_t sa, uint32_t sb){
                return sample_list[sa][major_axis] < sample_list[sb][major_axis];
            });

            real v_prev = -1e100;
            FOR_EACH(v, vertices) {
                dvec3 vv = data.vertices[v];
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
        FOR(a, vertices.size()) FOR(b, a) FOR(c, b) {
            plane candidate(data.vertices[vertices[a]], data.vertices[vertices[b]], data.vertices[vertices[c]]);
            evaluate_candidate(candidate, data, mesh, samples_begin, samples_end,
                    best_heuristic, best_hist, best_candidate);
        }
    }

    auto split = samples_begin;
    if (best_heuristic == std::numeric_limits<int64_t>::max()) {
        if (mesh.size() > 250)
            std::cout << "out of samples: " << mesh.size() << std::endl;
        // Out of samples -> split using new heuristic ignoring samples: abs(positive - negative) + 8 * splits
        FOR(a, vertices.size()) FOR(b, a) FOR(c, b) {
            plane candidate(data.vertices[vertices[a]], data.vertices[vertices[b]], data.vertices[vertices[c]]);

            Hist hist{ 0, 0, 0, 0 };
            FOR(f, mesh.size()) {
                const triangle3& face = data.faces[mesh[f]];

                constexpr real epsilon = 1e-5;
                int ai = classify(face.a, candidate, epsilon);
                int bi = classify(face.b, candidate, epsilon);
                int ci = classify(face.c, candidate, epsilon);

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

            int64_t heuristic = abs(hist.positive - hist.negative) + 8 * hist.stradle;
            if (heuristic < best_heuristic) {
                best_heuristic = heuristic;
                hist.overlap = mesh.size() - hist.positive - hist.negative - hist.stradle;
                best_hist = hist;
                best_candidate = candidate;
            }
        }
        split = samples_end = samples_begin;
    } else {
        split = std::partition(samples_begin, samples_end, [&data, &best_candidate](uint32_t s) {
            return best_candidate.distance(data.samples[s]) > 0;
        });
    }

    const auto psamples = std::distance(samples_begin, split);
    const auto nsamples = std::distance(split, samples_end);

    std::vector<uint32_t> pmesh, nmesh;
    pmesh.reserve(best_hist.positive + best_hist.stradle);
    nmesh.reserve(best_hist.negative + best_hist.stradle);
    FOR_EACH(f, mesh) {
        const triangle3& face = data.faces[f];

        constexpr real epsilon = 1e-5;
        int ai = classify(face.a, best_candidate, epsilon);
        int bi = classify(face.b, best_candidate, epsilon);
        int ci = classify(face.c, best_candidate, epsilon);

        if (ai > 0 || bi > 0 || ci > 0)
            pmesh.push_back(f);
        if (ai < 0 || bi < 0 || ci < 0)
            nmesh.push_back(f);
    }

    dvec3 vmin = data.vertices[vertices[0]];
    dvec3 vmax = data.vertices[vertices[0]];
    FOR_EACH(iv, vertices) {
        vmin = min(vmin, data.vertices[iv]);
        vmax = max(vmax, data.vertices[iv]);
    }
    dvec3 bs = vmax - vmin;

    uint32_t i = add_node();
    release(vertices);

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
    real score = 1 + (psamples * presult.second + nsamples * nresult.second) / samples;
    box_size[i] = bs;
    return std::pair<uint32_t, real>(i, score);
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

void SolidBSPTree::print(uint32_t n, int depth) {
    FOR(i, depth)
        std::cout << "  ";
    const uint32_t Inside = 0xFFFFFFFF;
    const uint32_t Outside = 0;

    if (n == Outside && depth != 0) {
        std::cout << "outside" << std::endl;
        return;
    }
    if (n == Inside) {
        std::cout << "inside" << std::endl;
        return;
    }
    std::cout << tfm::format("plane %.3f%% [%d %d %d %d] (%s)", 100 * percent[n],
        (int)hist[n].positive, (int)hist[n].negative, (int)hist[n].overlap, (int)hist[n].stradle, box_size[n]);
    std::cout << std::endl;
    print(node[n].positive, depth + 1);
    print(node[n].negative, depth + 1);
}

bool intersects_convex(const dvec3& v, const Mesh3d& mesh) {
    FOR_EACH(face, mesh)
        if (plane::sign(face.a, face.b, face.c, v) > 0)
            return false;
    return true;
}

bool intersects(const dvec3& v, const triangle3* mesh_begin, const triangle3* mesh_end) {
	// Can't just take the sdist sign of min_dist as two faces that share the edge can have the same dist but different sgn(sdist)
	real min_dist = std::numeric_limits<real>::max(), max_sdist = 0;
	constexpr real eps = 1e-5; // ??? arbitrary, need to account for rounding errors when computing distance(vertex, face)
    for (auto m = mesh_begin; m != mesh_end; m += 1) {
        const triangle3& face = *m;
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
	return max_sdist <= 0;
}

bool intersects(const dvec3& v, const Mesh3d& mesh) {
    return ::intersects(v, &mesh[0], &mesh[0] + mesh.size());
}

// TODO move to primitives.hh
template<typename RND>
dvec3 box_uniform_random(dvec3 vmin, dvec3 vmax, RND& rnd) {
    dvec3 d = vmax - vmin;
    real dmax = max(d.x, d.y, d.z);
    while (true) {
        dvec3 v = random_vector(rnd, 0, 1) * dmax;
        if (v.x <= d.x && v.y <= d.y && v.z <= d.z)
            return vmin + v;
    }
}

bool SolidBSPTree::intersects(const dvec3& v) {
    const uint32_t Inside = 0xFFFFFFFF;
    const uint32_t Outside = 0;

    uint32_t i = 0;
    while (true) {
        const Node& n = node[i];
        i = (n.divider.distance(v) > 0) ? n.positive : n.negative;
        if (i == Outside)
			return false;
		if (i == Inside)
			return true;
	}
}

SolidBSPTree::SolidBSPTree(const Mesh3d& mesh, uint32_t num_samples, std::default_random_engine& rnd) {
	BuildData data;
    data.rnd = rnd;

    dvec3 vmin = min(mesh);
    dvec3 vmax = max(mesh);

	// Init samples
	data.samples.resize(num_samples);
	std::vector<uint32_t> isamples(num_samples);
	FOR(i, num_samples) {
        data.samples[i] = box_uniform_random(vmin, vmax, data.rnd);
        isamples[i] = i;
	}

	// Init vertices and imesh
	std::unordered_map<dvec3, uint16_t> vec_index;
	std::vector<uint32_t> imesh;
	FOR_EACH(face, mesh) {
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
    std::cout << "TreeScore " << result.second << std::endl;
    print(0);

    int ab = 0, aa = 0, bb = 0;
    real n = num_samples;
    std::vector<dvec3> tester;
    tester.resize(n);
    FOR(j, n)
        tester[j] = box_uniform_random(vmin, vmax, data.rnd);

    FOR(j, n) {
        dvec3 v = tester[j];
        bool insideA = ::intersects(v, mesh);
        bool insideB = intersects(v);
        if (insideA)
            aa += 1;
        if (insideB)
            bb += 1;
        if (insideA == insideB)
            ab += 1;
    }
    dvec3 d = vmax - vmin;
    std::cout << tfm::format("Same %.05f, A %.05f, Tree %.05f Volume %.05f", ab / n, aa / n, bb / n, volume(mesh) / (d.x * d.y * d.z)) << std::endl;
}*/
