#include <core/diff.h>

static void DFS(PDiff a, unordered_set<Diff*>& visited, vector<PDiff>& out) {
    if (visited.count(a.get())) return;
    visited.insert(a.get());
    for (PDiff b : a->Inputs())
        if (b) DFS(b, visited, out);
    out.push_back(a);
}

vector<PDiff> TopoSort(const cspan<PDiff> heads) {
    vector<PDiff> out;
    unordered_set<Diff*> visited;
    for (PDiff a : heads) DFS(a, visited, out);
    return out;
}

bool IsBroadcastable(tensor_shape a, tensor_shape b) {
    return a.size != 0 && (a.volume == 1 || a == b.last(a.size));
}

struct BroadcastS : public Diff1 {
    BroadcastS(PDiff a, tensor_shape b) : Diff1(a) { Reshape(b); }
    void Forward() override { EACH(v) v[i] = va[0]; }
    void Backward() override { if (ga) EACH(g) ga[0] += g[i]; }
};

struct BroadcastT : public Diff1 {
    BroadcastT(PDiff a, tensor_shape b) : Diff1(a) { Reshape(b); }
    void Forward() override { EACH(v) v[i] = va[i % va.size]; }
    void Backward() override { if (ga) EACH(g) ga[i % ga.size] += g[i]; }
};

PDiff Broadcast(PDiff a, tensor_shape b) {
    if (a->shape == b) return a;
    if (a->size == 1) return make_shared<BroadcastS>(a, b);
    if (a->shape == b.last(a->rank)) return make_shared<BroadcastT>(a, b);
    Check(false);
    return nullptr;
}

MaxPool2D::MaxPool2D(PDiff a) : Diff1(a) {
    Check(a->rank == 2);
    const uint m = (a->shape(0) + 1) / 2;
    const uint n = (a->shape(1) + 1) / 2;
    Reshape({m, n});
}

void MaxPool2D::Forward() {
    // TODO edge condition on last row / column
    Check(a->shape(0) % 2 == 0);
    Check(a->shape(1) % 2 == 0);

    const uint m = shape(0);
    const uint n = shape(1);
    for (uint i = 0; i < m; i++) {
        for (uint j = 0; j < n; j++) {
            const uint p = i * 2, q = j * 2;
            v(i, j) = max(va(p, q), va(p + 1, q), va(p, q + 1), va(p + 1, q + 1));
        }
    }
}

void MaxPool2D::Backward() {
    if (!ga) return;
    Check(false);
    const uint m = shape(0);
    const uint n = shape(1);
    for (uint i = 0; i < m; i++) {
        for (uint j = 0; j < n; j++) {
            const uint p = i * 2, q = j * 2;
            g(i, j) = max(va(p, q), va(p + 1, q), va(p, q + 1), va(p + 1, q + 1));
        }
    }
}

void Forward(span<PDiff> nodes) {
    for (PDiff p : nodes) {
        Timestamp begin;
        p->Forward();
        Timestamp end;
        p->forward_ticks += begin.elapsed(end);
    }
}

void ResetGradients(span<PDiff> nodes) {
    for (PDiff p : nodes) EACH(p->g) p->g[i] = 0;
}

void Backward(span<PDiff> nodes) {
    for (auto i : range(nodes.size())) {
        auto p = nodes[nodes.size() - 1 - i];
        Timestamp begin;
        p->Backward();
        Timestamp end;
        p->backward_ticks += begin.elapsed(end);
    }
}

void GradientDescent(span<PDiff> nodes, const float alpha) {
    for (PDiff p : nodes)
        if (IsParam(p)) {
            Timestamp begin;
            EACH(p->v) p->v[i] -= alpha * p->g[i];
            Timestamp end;
            p->descend_ticks += begin.elapsed(end);
        }
}

bool Bounded(cspan<PDiff> nodes, const tensor::type limit) {
    for (PDiff p : nodes) {
        for (auto e : p->v) if (!std::isfinite(e) || abs(e) > limit) return false;
        for (auto e : p->g) if (!std::isfinite(e) || abs(e) > limit) return false;
    }
    return true;
}

auto ComputeIds(cspan<PDiff> nodes) {
    map<Diff*, string> ids;
    for (PDiff p : nodes) {
        if (IsConst(p) && p->v.size == 1 && p->name.empty()) {
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
    Aggregates<tensor::type> a(v.data(), v.data() + v.size);
    return format("(%s %s %s) %s", a.min, a.mean, a.max, sqrt(a.variance));
}

void Print(cspan<PDiff> nodes) {
    auto ids = ComputeIds(nodes);

    ulong ftotal = 0, btotal = 0, dtotal = 0;
    for (PDiff p : nodes)
        if (!(IsConst(p) && p->v.size == 1 && p->name.empty())) {
            ftotal += p->forward_ticks;
            btotal += p->backward_ticks;
            dtotal += p->descend_ticks;
        }
    ulong total = ftotal + btotal + dtotal;

    vector<string> table = {"id|type|inputs|shape|forward|backward|descend|ratio|values|gradients|"};
    string os;
    for (PDiff p : nodes)
        if (!(IsConst(p) && p->v.size == 1 && p->name.empty())) {
            os.clear();

            // id
            format_s(os, "%s|", ids.at(p.get()));

            // type
            string type = TypeName(*p.get());
            if (IsParam(p)) type = "Param";
            if (IsConst(p)) type = "Const";
            if (type.back() == 'T') type.pop_back();
            format_s(os, "%s|", type);

            // inputs
            for (PDiff e : p->Inputs())
                if (e) format_s(os, "%s ", ids.at(e.get()));
            os += '|';

            // shape
            format_s(os, "%s|", string(p->shape));

            // forward
            format_s(os, "%s|", (p->forward_ticks + 500) / 1000);

            // backward
            format_s(os, "%s|", (p->backward_ticks + 500) / 1000);

            // descend
            format_s(os, "%s|", (p->descend_ticks + 500) / 1000);

            // ratio
            ulong ticks = p->forward_ticks + p->backward_ticks + p->descend_ticks;
            format_s(os, "%.3f|", 100.0f * ticks / total);

            // values
            size_t summary = 6;
            format_s(os, "%s|", (p->v.size >= summary) ? Summary(p->v) : string(p->v));

            // gradients
            format_s(os, "%s|", (p->g.size >= summary) ? Summary(p->g) : string(p->g));

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

    // descend
    format_s(os, "%s|", (dtotal + 500) / 1000);

    // ratio
    format_s(os, "%.3f|", 100.0f);

    // values
    os += '|';

    // gradients
    os += '|';

    table << os;

    PrintTable(table, '|', " ", {4, 5, 6, 7});
}

Model::Model(PDiff loss, PDiff accuracy) : loss(loss), accuracy(accuracy) {
    Check(loss->size == 1);
    Check(accuracy->size == 1);
    // Recompute TopoSort as Setup can add more nodes.
    nodes = TopoSort({loss, accuracy});
    // TODO optimize it here:
    // - remove redundant diffs
    // - replace more expensive diffs with cheaper diffs
    // - fuse diffs
    // - if fanout of diff is 1 then its downstream diff can use = instead of += for gradients
    // - don't compute gradients which are not used!
    // - cpu vectorized kernels
    // - gpu vectorized kernels
    // - cpu multi-core kernels
}

vector<uint32_t> ShuffledInts(uint32_t size, std::mt19937_64& random) {
    vector<uint32_t> out(size);
    for (auto i : range(size)) out[i] = i;
    std::shuffle(out.begin(), out.end(), random);
    return out;
}

template <typename Func>
ulong Duration(const Func& func) {
    Timestamp begin;
    func();
    Timestamp end;
    return begin.elapsed(end);
}

Metrics Model::Epoch(cspan<pair<PDiff, tensor>> data, const float alpha, std::mt19937_64& random, bool verbose, uint epoch) {
    Check(data.size() > 0);
    const auto B = data[0].first->shape(0);
    const auto N = data[0].second.shape()[0];
    Check(N % B == 0, format("N %% B must be 0. N:%s B:%s", N, B));

    for (const auto& [key, value] : data) {
        Check(key->shape(0) == B);
        Check(value.shape()[0] == N);
        Check(value.shape().pop_front() == key->shape->pop_front());
    }

    auto samples = ShuffledInts(N, random);

    string msg;
    Accumulator<float> a_loss, a_accuracy, a_f_ticks, a_b_ticks, a_gd_ticks;
    ulong f_ticks = 0, b_ticks = 0, gd_ticks = 0;
    ulong message_ticks = 0;
    for (size_t i = 0; i < N; i += B) {
        for (const auto& [key, value] : data) {
            for (size_t j = 0; j < B; j++) key->v.slice(j).copy_from(value.slice(samples[i + j]));
        }

        f_ticks += Duration([&]() { Forward(); });

        a_loss << loss->v[0];
        a_accuracy << accuracy->v[0];

        if (alpha != 0) {
            ResetGradients();
            loss->g[0] = 1;
            b_ticks += Duration([&]() { Backward(); });
            gd_ticks += Duration([&]() { GradientDescent(alpha); });
        }

        a_f_ticks << f_ticks;
        a_b_ticks << b_ticks;
        a_gd_ticks << gd_ticks;

        ulong ticks = f_ticks + b_ticks + gd_ticks;
        if (ticks - message_ticks > long(1e10) && verbose) {
            message_ticks = ticks;

            for (char& c : msg) c = '\r';
            cout << msg;
            msg.clear();

            format_s(msg, "%s: %s/%s", epoch, i + B, N);
            format_s(msg, " loss:%.4f", a_loss.mean());
            format_s(msg, " acc:%.4f", a_accuracy.mean());
            msg += "   ";
            cout << msg;
            cout.flush();
        }
    }

    if (verbose) {
        for (char& c : msg) c = '\r';
        cout << msg;
        msg.clear();

        print("%s: %s/%s", epoch, N, N);
        print(" loss:%.4f", a_loss.mean());
        print(" acc:%.4f", a_accuracy.mean());
        println(" f_ticks:%h b_ticks:%h gd_ticks:%h", ulong(a_f_ticks.mean()), ulong(a_b_ticks.mean()),
                ulong(a_gd_ticks.mean()));
    }
    return {{"loss", a_loss.mean()}, {"accuracy", a_accuracy.mean()}};
}
