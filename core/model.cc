#include <core/model.h>

static void DFS(Diff a, unordered_set<DiffT*>& visited, vector<Diff>& out) {
    if (a == nullptr || visited.count(a.get())) return;
    visited.insert(a.get());
    for (Diff b : a->Inputs())
        if (b) DFS(b, visited, out);
    out.push_back(a);
}

vector<Diff> TopoSort(const cspan<Diff> heads) {
    vector<Diff> out;
    unordered_set<DiffT*> visited;
    for (Diff a : heads) DFS(a, visited, out);
    return out;
}

void Model::BeginEpoch(bool training) {
    for (auto p : m_begin_epoch_nodes) p->BeginEpoch(training);
}

void Model::Forward(bool training) {
    if (training) {
        for (auto p : m_forward_nodes) {
            Timestamp ts;
            p->Forward(true);
            p->forward_ticks += ts.elapsed();
        }
        return;
    }

    for (auto p : m_forward_nodes) if (p->inference) {
        Timestamp ts;
        p->Forward(false);
        p->forward_ticks += ts.elapsed();
    }
}

void Model::Backward(Diff loss) {
    for (auto p : m_params) {
        Timestamp ts;
        EACH(p->v) p->g[i] = 0;
        p->backward_ticks += ts.elapsed();
    }
    for (auto p : m_backward_nodes) {
        Timestamp ts;
        EACH(p->v) p->g[i] = 0;
        p->backward_ticks += ts.elapsed();
    }

    loss->g[0] = 1;

    for (auto p : m_backward_nodes) {
        Timestamp ts;
        p->Backward();
        p->backward_ticks += ts.elapsed();
    }

    optimizer->Optimize(span<ParamT*>(m_params));
}

void Model::EndEpoch(bool training) {
    for (auto p : m_end_epoch_nodes) p->EndEpoch(training);
}

bool Bounded(cspan<Diff> nodes, const tensor::type limit) {
    for (Diff p : nodes) {
        for (auto e : p->v)
            if (!std::isfinite(e) || abs(e) > limit) return false;
        for (auto e : p->g)
            if (!std::isfinite(e) || abs(e) > limit) return false;
    }
    return true;
}

auto ComputeIds(cspan<Diff> nodes) {
    map<DiffT*, string> ids;
    for (Diff p : nodes) {
        if (IsConst(p) && p->v.elements() == 1 && p->name.empty()) {
            ids.emplace(p.get(), format("%s", p->v[0]));
            continue;
        }

        string id = p->name.empty() ? "#0" : p->name;
        if (contains_value(ids, id)) {
            if (id == "#0") id = "#";
            int c = 1;
            while (contains_value(ids, format("%s%s", id, c))) c++;
            id = format("%s%s", id, c);
        }
        ids.emplace(p.get(), id);
    }
    return ids;
}

string Summary(const tensor v) {
    Aggregates<tensor::type> a(v.data(), v.data() + v.elements());
    return format("(%s %s %s) %s", a.min, a.mean, a.max, sqrt(a.variance));
}

// TODO show number of gradients in g with very tiny absolute values
void Model::Print() const {
    auto ids = ComputeIds(m_nodes);

    ulong ftotal = 0, btotal = 0, ototal = 0, rtotal = 0;
    for (Diff p : m_nodes)
        if (!(IsConst(p) && p->v.elements() == 1 && p->name.empty())) {
            ftotal += p->forward_ticks;
            btotal += p->backward_ticks;
        }
    const ulong total = ftotal + btotal;

    vector<string> table = {"id|type|inputs|shape|forward|backward|ratio|values|gradients|"};
    string os;
    for (Diff p : m_nodes)
        if (!(IsConst(p) && p->v.elements() == 1 && p->name.empty())) {
            os.clear();

            // id
            format_s(os, "%s|", ids.at(p.get()));

            // type
            string type = TypeName(*p.get());
            if (IsConst(p)) type = "Const";
            if (type.back() == 'T') type.pop_back();
            format_s(os, "%s|", type);

            // inputs
            for (Diff e : p->Inputs())
                if (e) format_s(os, "%s ", ids.at(e.get()));
            os += '|';

            // shape
            format_s(os, "%s|", string(p->shape()));

            // forward
            format_s(os, "%s|", (p->forward_ticks + 500) / 1000);

            // backward
            format_s(os, "%s|", (p->backward_ticks + 500) / 1000);

            // ratio
            ulong ticks = p->forward_ticks + p->backward_ticks;
            format_s(os, "%.3f|", 100.0f * ticks / total);

            // values
            size_t summary = 9;
            format_s(os, "%s|", (p->v.elements() > summary) ? Summary(p->v) : string(p->v));

            // gradients
            format_s(os, "%s|", (p->g.elements() > summary) ? Summary(p->g) : string(p->g));

            table << os;
        }

    os.clear();

    // id
    os += '|';

    // type
    os += '|';

    // inputs
    os += '|';

    // shape
    os += '|';

    // forward
    format_s(os, "%s|", (ftotal + 500) / 1000);

    // backward
    format_s(os, "%s|", (btotal + 500) / 1000);

    // ratio
    format_s(os, "%.3f|", 100.0f);

    // values
    os += '|';

    // gradients
    os += '|';

    table << os;

    std::unique_lock lock(g_cout_mutex);
    PrintTable(table, '|', " ", {4, 5, 6});
}

void Model::Iterate(size_t iterations, Diff loss) {
    for (size_t i : range(iterations)) {
        BeginEpoch(true);
        Forward(true);
        Backward(loss);
        EndEpoch(true);
    }
}

bool HasBackward(Diff p) {
    DiffT::has_overload = true;
    p->Backward();
    return DiffT::has_overload;
}

Model::Model(std::initializer_list<Diff> heads) : m_nodes(TopoSort(heads)) {
    // TODO optimize it here:
    // - remove redundant diffs
    // - replace more expensive diffs with cheaper diffs
    // - fuse diffs
    // - if fanout of diff is 1 then its downstream diff can use = instead of += for gradients
    // - don't compute gradients which are not used!
    // - cpu vectorized kernels
    // - gpu vectorized kernels
    // - cpu multi-core kernels

    for (Diff p : m_nodes) {
        DiffT::has_overload = true;
        p->Forward(false);
        if (DiffT::has_overload) m_forward_nodes.push_back(p.get());
    }

    for (Diff p : m_nodes) {
        if (HasBackward(p)) m_backward_nodes.push_back(p.get());
    }
    std::reverse(m_backward_nodes.begin(), m_backward_nodes.end());

    for (Diff p : m_nodes) {
        DiffT::has_overload = true;
        p->BeginEpoch(false);
        if (DiffT::has_overload) m_begin_epoch_nodes.push_back(p.get());
    }

    for (Diff p : m_nodes) {
        DiffT::has_overload = true;
        p->EndEpoch(false);
        if (DiffT::has_overload) m_end_epoch_nodes.push_back(p.get());
    }

    for (Diff p : m_nodes) {
        auto param = dynamic_cast<ParamT*>(p.get());
        if (param) m_params.push_back(param);
    }
}

Metrics Model::Epoch(Diff loss, Diff accuracy, cspan<pair<Diff, tensor>> data, std::mt19937_64& random, bool verbose, uint epoch, bool backward) {
    Check(data.size() > 0);
    const dim_t B = data[0].first->dim(0);
    const dim_t N = data[0].second.dim(0);
    Check(N % B == 0, format("N %% B must be 0. N:%s B:%s", N, B));

    for (const auto& [key, value] : data) {
        Check(key->dim(0) == B);
        Check(value.dim(0) == N);
        Check(value.shape().pop_front() == key->shape().pop_front(), format("%s vs %s", string(value.shape()), string(key->shape())));
    }

    if (m_samples.size() != N) {
        m_samples.resize(N);
        for (auto i : range(N)) m_samples[i] = i;
    }
    std::shuffle(m_samples.begin(), m_samples.end(), random);

    string msg;
    Accumulator<float> a_loss, a_accuracy;
    Accumulator<float> a_f_ticks, a_b_ticks;
    ulong f_ticks = 0, b_ticks = 0;
    ulong message_ticks = 0;

    for (DiffT* p : m_begin_epoch_nodes) p->BeginEpoch(true);

    for (size_t i = 0; i < N; i += B) {
        for (const auto& [key, value] : data) {
            FOR(j, B) key->v.slice(j).copy_from(value.slice(m_samples[i + j]));
        }

        f_ticks += Duration([&]() { Forward(true); });
        if (backward) b_ticks += Duration([&]() { Backward(loss); });
        a_f_ticks << f_ticks;
        a_b_ticks << b_ticks;

        a_loss << loss->v[0];
        if (accuracy != nullptr) a_accuracy << accuracy->v[0];

        ulong ticks = f_ticks + b_ticks;
        if (ticks - message_ticks > long(1e10) && verbose) {
            message_ticks = ticks;

            for (char& c : msg) c = '\r';
            cout << msg;
            msg.clear();

            format_s(msg, "%s: %s/%s", epoch, i + B, N);
            format_s(msg, " loss:%.5f", a_loss.mean());
            if (accuracy != nullptr) format_s(msg, " accuracy:%.5f", a_accuracy.mean());
            msg += "   ";
            cout << msg;
            cout.flush();
        }
    }

    for (DiffT* p : m_end_epoch_nodes) p->EndEpoch(true);

    if (verbose) {
        for (char& c : msg) c = '\r';
        cout << msg;
        msg.clear();

        print("%s: %s/%s", epoch, N, N);
        print(" loss:%.5f", a_loss.mean());
        if (accuracy) print(" accuracy:%.5f", a_accuracy.mean());
        print(" f_ticks:%h", ulong(a_f_ticks.mean()));
        print(" b_ticks:%h", ulong(a_b_ticks.mean()));
        println();
    }

    Metrics metrics = {{"loss", a_loss.mean()}};
    if (accuracy) metrics.emplace("accuracy", a_accuracy.mean());
    return metrics;
}

void NormalizeDataset(tensor a) {
    vtensor q(a.shape().pop_front());

    const auto N = a.shape()[0];
    for (auto s : range(N)) {
        auto sample = a.slice(s);
        EACH(q) q[i] += sample[i];
    }
    EACH(q) q[i] /= N;
    for (auto s : range(N)) {
        auto sample = a.slice(s);
        EACH(q) sample[i] -= q[i];
    }

    EACH(q) q[i] = 0;
    for (auto s : range(N)) {
        auto sample = a.slice(s);
        EACH(q) q[i] += sqr(sample[i]);
    }
    EACH(q) q[i] = 1 / sqrt(q[i] / N + 1e-8);
    for (auto s : range(N)) {
        auto sample = a.slice(s);
        EACH(q) sample[i] *= q[i];
    }
}
