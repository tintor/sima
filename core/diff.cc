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

    string os = "id|type|inputs|shape|";
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

            if (values) format_s(os, "%s|", string(p->v));

            if (gradients) format_s(os, "%s|", string(p->g));

            table << os;
        }
    PrintTable(table, '|', " ");
}

void Model::Print() const {
    Check(compiled);
    ::Print(nodes);
}

vector<uint32_t> ShuffledInts(uint32_t size, std::mt19937_64& random) {
    vector<uint32_t> out(size);
    for (auto i : range(size)) out[i] = i;
    std::shuffle(out.begin(), out.end(), random);
    return out;
}

Metrics Model::Epoch(const tensor in, const tensor ref, const optional<float> alpha) {
    Check(compiled);
    Check(in.shape().size() > 0);
    Check(ref.shape().size() > 0);
    Check(in.shape()[0] == ref.shape()[0]);

    Check(IsData(input));
    Check(IsData(reference));
    Check(in.shape().pop_front() == input->shape().pop_front());
    Check(ref.shape().pop_front() == reference->shape().pop_front());

    std::random_device rd;
    std::mt19937_64 random(rd());

    const size_t n = in.shape()[0];
    auto samples = ShuffledInts(n, random);

    string msg;
    Accumulator<float> mean_loss, mean_accuracy;
    for (auto i : range(n)) {
        auto index = samples[i];
        Copy(in.slice(index), input->v);
        Copy(ref.slice(index), reference->v);
        Forward();

        mean_loss << loss->v[0];
        mean_accuracy << accuracy->v[0];

        if (alpha) {
            ResetGradients();
            loss->g[0] = 1;
            Backward();
            GradientDescent(*alpha);
        }

        for (char& c : msg) c = '\r';
        cout << msg;
        msg.clear();

        format_s(msg, "Epoch: %f%%", 100. * (i + 1) / n);
        format_s(msg, " %s/%s", i + 1, n);
        format_s(msg, " loss:%.4f", mean_loss.mean());
        format_s(msg, " acc:%.4f", mean_accuracy.mean());
        cout << msg;
        cout.flush();
    }
    print("\n");
    return {{"loss", mean_loss.mean()}, {"accuracy", mean_accuracy.mean()}};
}
