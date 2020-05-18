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

#define va (a->v)
#define vb (b->v)
#define vc (c->v)

#define ga (a->g)
#define gb (b->g)
#define gc (c->g)

bool IsBroadcastable(tensor_shape a, tensor_shape b) {
    return a.size() != 0 && (a.volume() == 1 || a == b.last(a.size()));
}

struct IBroadcastT : public DiffA {
    IBroadcastT(PDiff a, tensor_shape b) : DiffA(a) {
        Check(IsBroadcastable(a->shape(), b), format("%s %s", a->shape(), b));
        Reshape(b);
    }
    void Setup() override { Fail("forbidden"); }
    void Forward() override { EACH(v) v[i] = va[i % va.size()]; }
    void Backward() override { EACH(g) ga[i % ga.size()] += g[i]; }
};

PDiff IBroadcast(PDiff a, tensor_shape b) {
    return make_shared<IBroadcastT>(a, b);
}

void MaxPool2D::Setup() {
    Check(a->shape().size() == 2);
    // TODO overflow check
    const uint16_t m = (a->shape()[0] + 1) / 2;
    const uint16_t n = (a->shape()[1] + 1) / 2;
    Reshape({m, n});
}

void MaxPool2D::Forward() {
    // TODO edge condition on last row / column
    Check(a->shape()[0] % 2 == 0);
    Check(a->shape()[1] % 2 == 0);

    const size_t m = shape()[0];
    const size_t n = shape()[1];
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            const size_t p = i * 2, q = j * 2;
            v(i, j) = max(va(p, q), va(p + 1, q), va(p, q + 1), va(p + 1, q + 1));
        }
    }
}

void MaxPool2D::Backward() {
    if (!ga) return;
    Check(false);
    const size_t m = shape()[0];
    const size_t n = shape()[1];
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            const size_t p = i * 2, q = j * 2;
            g(i, j) = max(va(p, q), va(p + 1, q), va(p, q + 1), va(p + 1, q + 1));
        }
    }
}

#undef va
#undef vb
#undef vc

#undef ga
#undef gb
#undef gc

void Setup(span<PDiff> nodes) {
    for (PDiff p : nodes) p->Setup();
}

void InitBatches(span<PDiff> nodes, int batch_size) {
    for (PDiff p : nodes)
        if (IsData(p)) p->v.reshape(p->v.shape().set(0, batch_size));
}

void ResetTicks(span<PDiff> nodes) {
    for (PDiff p : nodes) p->forward_ticks = p->backward_ticks = 0;
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
    EACH(nodes) {
        auto p = nodes[nodes.size() - 1 - i];
        Timestamp begin;
        p->Backward();
        Timestamp end;
        p->backward_ticks += begin.elapsed(end);
    }
}

void GradientDescent(span<PDiff> nodes, const float alpha) {
    for (PDiff p : nodes)
        if (IsParam(p)) EACH(p->v) p->v[i] -= alpha * p->g[i];
}

void CheckBounded(cspan<PDiff> nodes, const tensor::type limit) {
    for (PDiff p : nodes) {
        EACH(p->v) Check(std::isfinite(p->v[i]) && abs(p->v[i]) <= limit, "value not bounded");
        EACH(p->g) Check(std::isfinite(p->g[i]) && abs(p->g[i]) <= limit, "gradient not bounded");
    }
}

auto ComputeIds(cspan<PDiff> nodes) {
    map<Diff*, string> ids;
    for (PDiff p : nodes) {
        if (IsConst(p) && p->v.size() == 1 && p->name.empty()) {
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
    Accumulator<tensor::type> acc;
    for (auto i : range(v.size())) acc << v[i];
    return format("(%s %s %s) %s", acc.min(), acc.mean(), acc.max(), acc.stdev());
}

void Print(cspan<PDiff> nodes) {
    auto ids = ComputeIds(nodes);
    vector<string> table = {"id|type|inputs|shape|fticks|bticks|values|gradients|"};
    string os;
    for (PDiff p : nodes)
        if (!(IsConst(p) && p->v.size() == 1 && p->name.empty())) {
            os.clear();

            // id
            format_s(os, "%s|", ids.at(p.get()));

            // type
            string type = TypeName(*p.get());
            if (IsParam(p)) type = "Param";
            if (IsConst(p)) type = "Const";
            if (IsData(p)) type = "Data";
            if (type.back() == 'T') type.pop_back();
            format_s(os, "%s|", type);

            // inputs
            for (PDiff e : p->Inputs())
                if (e) format_s(os, "%s ", ids.at(e.get()));
            os += '|';

            // shape
            format_s(os, "%s|", string(p->shape()));

            // fticks
            format_s(os, "%s|", (p->forward_ticks + 500) / 1000);

            // bticks
            format_s(os, "%s|", (p->backward_ticks + 500) / 1000);

            size_t summary = 6;
            format_s(os, "%s|", (p->v.size() >= summary) ? Summary(p->v) : string(p->v));

            format_s(os, "%s|", (p->g.size() >= summary) ? Summary(p->g) : string(p->g));

            table << os;
        }
    PrintTable(table, '|', " ");
}

Model::Model(PDiff loss, PDiff accuracy, int batch_size)
    : loss(loss), accuracy(accuracy), nodes(TopoSort({loss, accuracy})) {
    InitBatches(nodes, batch_size);
    Setup(nodes);
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

Metrics Model::Epoch(cspan<pair<PDiff, tensor>> data, const float alpha, const size_t seed) {
    Check(data.size() > 0);
    const auto B = data[0].first->shape()[0];
    const auto N = data[0].second.shape()[0];
    Check(N % B == 0);

    for (const auto& [key, value] : data) {
        Check(key->shape()[0] == B);
        Check(value.shape()[0] == N);
        Check(value.shape().pop_front() == key->shape().pop_front(), format("%s %s", value.shape(), key->shape()));
    }

    std::mt19937_64 random(seed);
    auto samples = ShuffledInts(N, random);

    string msg;
    Accumulator<float> a_loss, a_accuracy, a_f_ticks, a_b_ticks, a_gd_ticks;
    ulong f_ticks = 0, b_ticks = 0, gd_ticks = 0;
    for (size_t i = 0; i < N; i += B) {
        for (const auto& [key, value] : data) {
            for (size_t j = 0; j < B; j++)
                key->v.slice(j).copy_from(value.slice(samples[i + j]));
        }

        f_ticks += Duration([&](){ Forward(); });

        a_loss << loss->v[0];
        a_accuracy << accuracy->v[0];

        if (alpha != 0) {
            ResetGradients();
            loss->g[0] = 1;
            b_ticks += Duration([&](){ Backward(); });
            gd_ticks += Duration([&](){ GradientDescent(alpha); });
        }

        for (char& c : msg) c = '\r';
        cout << msg;
        msg.clear();

        a_f_ticks << f_ticks;
        a_b_ticks << b_ticks;
        a_gd_ticks << gd_ticks;

        format_s(msg, "%s/%s", i + B, N);
        //format_s(msg, " loss:%.4f", a_loss.mean());
        //format_s(msg, " acc:%.4f", a_accuracy.mean());
        //msg += "   ";
        cout << msg;
        cout.flush();

        //::Print(nodes);
        // CheckBounded(nodes, 1e6);
    }

    for (char& c : msg) c = '\r';
    cout << msg;
    msg.clear();

    print("%s/%s", N, N);
    print(" loss:%.4f", a_loss.mean());
    print(" acc:%.4f", a_accuracy.mean());
    println(" f_ticks:%h b_ticks:%h gd_ticks:%h", ulong(a_f_ticks.mean()), ulong(a_b_ticks.mean()), ulong(a_gd_ticks.mean()));
    return {{"loss", a_loss.mean()}, {"accuracy", a_accuracy.mean()}};
}
