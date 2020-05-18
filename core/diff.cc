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

void Print(cspan<PDiff> nodes, bool values, bool gradients) {
    auto ids = ComputeIds(nodes);
    vector<string> table;

    string os = "id|type|inputs|shape|fticks|bticks|";
    if (values) os += "values|";
    if (gradients) os += "gradients|";
    table << os;

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

            if (values) format_s(os, "%s|", string(p->v));

            if (gradients) format_s(os, "%s|", string(p->g));

            table << os;
        }
    PrintTable(table, '|', " ");
}

vector<uint32_t> ShuffledInts(uint32_t size, std::mt19937_64& random) {
    vector<uint32_t> out(size);
    for (auto i : range(size)) out[i] = i;
    std::shuffle(out.begin(), out.end(), random);
    return out;
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
    Accumulator<float> mean_loss, mean_accuracy;
    for (size_t i = 0; i < N; i += B) {
        for (const auto& [key, value] : data) {
            for (size_t j = 0; j < B; j++)
                key->v.slice(j).copy_from(value.slice(samples[i + j]));
        }

        Forward();

        mean_loss << loss->v[0];
        mean_accuracy << accuracy->v[0];

        if (alpha != 0) {
            ResetGradients();
            loss->g[0] = 1;
            Backward();
            GradientDescent(alpha);
        }

        for (char& c : msg) c = '\r';
        cout << msg;
        msg.clear();

        format_s(msg, "%s/%s", i + 1, N);
        format_s(msg, " loss:%.4f", mean_loss.mean());
        format_s(msg, " acc:%.4f", mean_accuracy.mean());
        cout << msg;
        cout.flush();

        // CheckBounded(nodes, 1e6);
    }
    println();
    return {{"loss", mean_loss.mean()}, {"accuracy", mean_accuracy.mean()}};
}
