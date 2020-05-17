#include <core/diff.h>

#include <cxxabi.h>

static void DFS(PDiff a, unordered_set<Diff*>& visited, vector<PDiff>& out) {
    if (visited.count(a.get())) return;
    for (PDiff b : a->Inputs()) if (b) DFS(b, visited, out);
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

template<typename T>
inline string TypeName(const T& value) {
    int status;
    char* demangled = abi::__cxa_demangle(typeid(value).name(), 0, 0, &status);
    ON_SCOPE_EXIT(free(demangled));
    return demangled;
}

void Model::Print() const {
    Check(compiled);
    map<Diff*, string> ids;
    for (PDiff p : nodes) {
        if (IsConstant(p) && p->v.size() == 1 && p->name.empty()) {
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

        print("%s = %s", id, TypeName(*p.get()));
        for (PDiff p : p->Inputs()) if (p) print(" %s", ids.at(p.get()));
        print(" %s\n", string(p->shape()));
    }
}

vector<uint32_t> ShuffledInts(uint32_t size, std::mt19937_64& random) {
    vector<uint32_t> out(size);
    for (auto i : range(size)) out[i] = i;
    std::shuffle(out.begin(), out.end(), random);
    return out;
}

Metrics SGD::Train(Model& model, const tensor in, const tensor ref) {
    Check(model.compiled);
    Check(in.shape().size() > 0);
    Check(ref.shape().size() > 0);
    Check(in.shape()[0] == ref.shape()[0]);

    Check(IsData(model.input));
    Check(IsData(model.reference));
    Check(in.shape().pop_front() == model.input->shape().pop_front());
    Check(ref.shape().pop_front() == model.reference->shape().pop_front());

    const size_t n = in.shape()[0];
    auto samples = ShuffledInts(n, random);

    string message;
    Accumulator<float> loss_accum, accuracy_accum;
    for (auto i : range(n)) {
        auto index = samples[i];
        Copy(in.slice(index), model.input->v);
        Copy(ref.slice(index), model.reference->v);
        model.Forward();

        auto loss = model.loss->v[0];
        auto accuracy = model.accuracy->v[0];
        loss_accum << loss;
        accuracy_accum << accuracy;

        model.Backward();
        model.Adjust(alpha);

        for (char& c : message) c = '\r';
        print("%s", message);
        message = format("Epoch: %f%% %s/%s loss:%.4f acc:%.4f",
            100. * (i + 1) / n, i + 1, n, loss_accum.mean(), accuracy_accum.mean());
        print("%s", message);
        cout.flush();
    }
    print("\n");
    return {};
}

Metrics SGD::Test(Model& model, const tensor in, const tensor ref) {
    return {};
}
