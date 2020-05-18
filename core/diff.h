#pragma once
#include <core/auto.h>
#include <core/range.h>
#include <core/std.h>
#include <core/tensor.h>

// TODO Data augmentation operators:
// - non-gradient batch norm
// - non-gradient 2d mirror
// - non-gradient 2d small random translate
// - non-gradient 2d small random rotate
// - non-gradient 2d small random scale

struct Diff;

using PDiff = std::shared_ptr<Diff>;

#define va (a->v)
#define vb (b->v)
#define vc (c->v)

#define ga (a->g)
#define gb (b->g)
#define gc (c->g)

struct Diff;

struct Diff {
    virtual array<PDiff, 3> Inputs() { return {nullptr, nullptr, nullptr}; }
    virtual void Setup() {}
    virtual void Forward() {}
    virtual void Backward() {}

    void Reshape(tensor_shape s) {
        v.reshape(s);
        g.reshape(s);
    }

    auto shape() const { return v.shape(); }
    size_t size() const { return v.size(); }

    string name;
    vtensor v, g;

    ulong forward_ticks = 0;
    ulong backward_ticks = 0;
};

#define EACH(V) for (auto i : range(V.size()))

#define Declare1(Func) \
    inline PDiff Func(PDiff a) { return make_shared<Func##T>(a); }
#define Declare2(Func, Type) \
    inline PDiff Func(PDiff a, PDiff b) { return make_shared<Type>(a, b); }

// Returns a list of all diffs that need to be computed (in order) before all goal diffs.
// Goal diffs are also ordered if one depends on any other.
vector<PDiff> TopoSort(const cspan<PDiff> heads);

inline vector<PDiff> LoadModel(string_view filename) { return {}; }

inline void SaveModel(span<PDiff> model, string_view filename) {}

struct DiffA : public Diff {
    DiffA(PDiff a) : a(a) {}
    array<PDiff, 3> Inputs() override { return {a, nullptr, nullptr}; }
    void Setup() override { Reshape(a->shape()); }

    PDiff a;
};

bool IsBroadcastable(tensor_shape a, tensor_shape b);
PDiff IBroadcast(PDiff a, tensor_shape b);

struct DiffAB : public Diff {
    DiffAB(PDiff a, PDiff b) : a(a), b(b) {}
    array<PDiff, 3> Inputs() override { return {a, b, nullptr}; }

    void Setup() override {
        if (a->shape() != b->shape()) {
            if (a->size() > b->size()) b = IBroadcast(b, a->shape());
            if (a->size() < b->size()) a = IBroadcast(a, b->shape());
            Check((a->size() == 1 && b->size() == 1) ||  a->shape() == b->shape(), format("%s %s", a->shape(), b->shape()));
        }
        Reshape(a->shape());
    }

    PDiff a, b;
};

struct DiffABC : public Diff {
    DiffABC(PDiff a, PDiff b, PDiff c) : a(a), b(b), c(c) {}
    array<PDiff, 3> Inputs() override { return {a, b, c}; }

    PDiff a, b, c;
};

// Ignores gradients and batches.
inline PDiff Const(tensor::type c, string_view name = "") {
    auto p = make_shared<Diff>();
    p->v.reshape({1});
    p->v[0] = c;
    p->name = name;
    return p;
}

inline PDiff Constant(const tensor c, string_view name = "") {
    auto p = make_shared<Diff>();
    p->v = c;
    p->name = name;
    return p;
}

struct Init {
    virtual tensor::type get() { return 0; }
};

struct NormalInit : public Init {
    NormalInit(tensor::type variance, size_t seed) : dis(0, variance), random(seed) {}
    NormalInit(tensor::type variance) : NormalInit(variance, std::random_device()()) {}
    tensor::type get() override { return dis(random); }

    std::normal_distribution<tensor::type> dis;
    std::mt19937_64 random;
};

struct UniformInit : public Init {
    UniformInit(tensor::type min, tensor::type max, size_t seed) : dis(min, max), random(seed) {}
    tensor::type get() override { return dis(random); }

    std::uniform_real_distribution<tensor::type> dis;
    std::mt19937_64 random;
};

// Accepts gradients, ignores batches. Learnable parameter, can be saved.
inline PDiff Param(tensor_shape shape, string_view name, shared_ptr<Init> init = make_shared<Init>()) {
    auto p = make_shared<Diff>();
    Check(shape.size() > 0, "Parameter shape.size() must be non-zero");
    p->Reshape(shape);
    EACH(p->v) p->v[i] = init->get();
    p->name = name;
    return p;
}

// Input and reference value. Not saved.
// Ignores gradients, requires batches.
inline PDiff Data(tensor_shape shape, string_view name = "") {
    auto p = make_shared<Diff>();
    p->name = name;
    p->v.reshape(shape);
    return p;
}

struct ReluT : public DiffA {
    ReluT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = (va[i] > 0) ? va[i] : 0; }
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) ? g[i] : 0; }
};

struct LeakyReluT : public DiffA {
    LeakyReluT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (va[i] / 10); }
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) ? g[i] : (g[i] / 10); }
};

struct SigmoidT : public DiffA {
    SigmoidT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = 1 / (1 + exp(-va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i] * (1 - v[i]); }
};

struct SqrT : public DiffA {
    SqrT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = va[i] * va[i]; }
    void Backward() override { EACH(ga) ga[i] += g[i] * 2 * va[i]; }
};

struct SqrtT : public DiffA {
    SqrtT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = sqrt(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / 2 / v[i]; }
};

struct ExpT : public DiffA {
    ExpT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = exp(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i]; }
};

struct LogT : public DiffA {
    LogT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = log(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / va[i]; }
};

struct AbsT : public DiffA {
    AbsT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = abs(va[i]); }
    void Backward() override {
        EACH(ga) {
            if (va[i] > 0) ga[i] += g[i];
            if (va[i] < 0) ga[i] -= g[i];
        }
    }
};

struct CosT : public DiffA {
    CosT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = cos(va[i]); }
    void Backward() override { EACH(ga) ga[i] -= sin(va[i]); }
};

struct SinT : public DiffA {
    SinT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = sin(va[i]); }
    void Backward() override { EACH(ga) ga[i] += cos(va[i]); }
};

Declare1(Relu);
Declare1(LeakyRelu);
Declare1(Sigmoid);
Declare1(Sqr);
Declare1(Sqrt);
Declare1(Exp);
Declare1(Log);
Declare1(Abs);
Declare1(Cos);
Declare1(Sin);

#if 0
// Y = A * b + c
// b and c are scalars
struct BroadcastFMA : public DiffABC {
    BroadcastFMA(PDiff a, PDiff b, PDiff c) : DiffABC(a, b, c) {}

    void Setup() override {
        Check(b.size() == 1);
        Check(c.size() == 1);
        y.reshape(a.shape());
        dy.reshape(a.shape());
    }

    void Forward() override {
        for (auto i : range(v)()) v[i] = va[i] * b[0] + c[0];
    }

    void Backward() override {
        for (auto i : range(v)()) ga[i] += g[i] * b[0];
        for (auto i : range(v)()) gb[0] += g[i] * va[i];
        if (c.dy) for (auto i : range(v)()) c.g[0] += g[i];
    }
};
#endif


// [1] -> [N] = [N 1]
// [4] -> [N 4]
struct BroadcastT : public DiffAB {
    BroadcastT(PDiff a, PDiff b) : DiffAB(a, b) {}

    void Setup() override {
        Check(IsBroadcastable(a->shape(), b->shape()), format("%s %s", a->shape(), b->shape()));
        Reshape(b->shape());
    }

    void Forward() override { EACH(v) v[i] = va[i % va.size()]; }
    void Backward() override { if (ga) EACH(g) ga[i % ga.size()] += g[i]; }
};

struct Add : public DiffAB {
    Add(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] + vb[i]; }

    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] += g[i];
    }
};

struct Sub : public DiffAB {
    Sub(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] - vb[i]; }

    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] += -g[i];
    }
};

struct Mul : public DiffAB {
    Mul(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] * vb[i]; }

    void Backward() override {
        EACH(ga) ga[i] += g[i] * vb[i];
        EACH(gb) gb[i] += g[i] * va[i];
    }
};

struct Div : public DiffAB {
    Div(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] / vb[i]; }

    void Backward() override {
        EACH(ga) ga[i] += g[i] / vb[i];
        EACH(gb) gb[i] -= g[i] / (va[i] * va[i]);
    }
};

Declare2(Broadcast, BroadcastT);
Declare2(operator+, Add);
Declare2(operator-, Sub);
Declare2(operator*, Mul);
Declare2(operator/, Div);

#define RELATION(NAME, OP)                                                   \
    struct NAME : public DiffAB {                                            \
        NAME(PDiff a, PDiff b) : DiffAB(a, b) {}                             \
        void Forward() override { EACH(v) v[i] = (va[i] OP vb[i]) ? 1 : 0; } \
    };                                                                       \
    Declare2(operator OP, NAME)

RELATION(Greater, >);
RELATION(Less_, <);
RELATION(GreaterOrEqual, >=);
RELATION(LessOrEqual, <=);

#undef RELATION

// TODO Broadcast
#define CONST_OVERLOAD(OP)                                                     \
    inline auto operator OP(tensor::type a, PDiff b) { return Broadcast(Const(a), b) OP b; } \
    inline auto operator OP(PDiff a, tensor::type b) { return a OP Broadcast(Const(b), a); }

CONST_OVERLOAD(+);
CONST_OVERLOAD(-);
CONST_OVERLOAD(*);
CONST_OVERLOAD(/);
CONST_OVERLOAD(>);
CONST_OVERLOAD(<);
CONST_OVERLOAD(>=);
CONST_OVERLOAD(<=);

#undef CONST_OVERLOAD

struct MinT : public DiffAB {
    MinT(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = min(va[i], vb[i]); }

    void Backward() override {
        EACH(ga) {
            if (va[i] < vb[i]) ga[i] += g[i];
            if (va[i] == vb[i]) ga[i] += g[i] / 2;
        }
        EACH(gb) {
            if (vb[i] < va[i]) gb[i] += g[i];
            if (vb[i] == va[i]) gb[i] += g[i] / 2;
        }
    }
};

struct MaxT : public DiffAB {
    MaxT(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = max(va[i], vb[i]); }

    void Backward() override {
        EACH(ga) {
            if (va[i] > vb[i]) ga[i] += g[i];
            if (va[i] == vb[i]) ga[i] += g[i] / 2;
        }
        EACH(gb) {
            if (vb[i] > va[i]) gb[i] += g[i];
            if (vb[i] == va[i]) gb[i] += g[i] / 2;
        }
    }
};

Declare2(Min, MinT);
Declare2(Max, MaxT);

struct Concat : public DiffAB {
    Concat(PDiff a, PDiff b) : DiffAB(a, b) {}

    void Setup() override {
        // TODO generalize for more dimensions
        Check(a->shape().size() == 1);
        Check(b->shape().size() == 1);
        Reshape(tensor_shape(a->size() + b->size()));
    }

    void Forward() override {
        EACH(v) v[i] = va[i];
        EACH(v) v[i + va.size()] = vb[i];
    }

    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] += g[i + va.size()];
    }
};

// TODO Splice tensor operator

struct Conv1D : public DiffAB {
    Conv1D(PDiff a, PDiff b, int offset) : DiffAB(a, b), offset(offset) {}

    void Setup() override {
        Check(a->shape().size() == 1);
        Check(b->shape().size() == 1);
        Reshape(a->shape());
    }

    void Forward() override {
        EACH(v) {
            tensor::type sum = 0;
            // TODO check
            size_t begin = (offset > i) ? offset - i : 0;
            size_t end = min(vb.size(), offset - i + va.size());
            for (size_t j = begin; j < end; j++) sum += va[i + j - offset] * vb[j];
            v[i] = sum;
        }
    }

    void Backward() override { Check(false, "not implemented"); }

    int offset;
};

struct Deconv1D : public DiffAB {
    Deconv1D(PDiff a, PDiff b) : DiffAB(a, b) {}
};

struct MaxPool1D : public DiffA {
    MaxPool1D(PDiff a) : DiffA(a) {}

    void Setup() override {
        Check(a->shape().size() == 1);
        const uint16_t m = (a->shape()[0] + 1) / 2;
        Reshape({m});
    }

    void Forward() override {
        const size_t m = a->shape()[0] / 2;
        for (size_t i = 0; i < m; i++) v[i] = max(va[i * 2], va[i * 2 + 1]);
        if (a->shape()[0] % 2 == 1) v[m] = va[m * 2];
    }

    void Backward() override {
        if (!ga) return;
        Check(false);
        const size_t m = shape()[0];
        for (size_t i = 0; i < m; i++) g[i] = max(va[i * 2], va[i * 2 + 1]);
    }
};

struct MaxPool2D : public DiffA {
    MaxPool2D(PDiff a) : DiffA(a) {}
    void Setup() override;
    void Forward() override;
    void Backward() override;
};

struct Reshape : public DiffA {
    Reshape(PDiff a, tensor_shape shape) : DiffA(a) { v.reshape(shape); }

    void Setup() override {
        tensor_shape s = v.shape();
        if (s[0] == 0) s = s.set(0, 1);  // TODO mini_batch size
        v.reshape(s);
        Check(a->size() == v.size());
    }

    void Forward() override { EACH(v) v[i] = va[i]; }
    void Backward() override { EACH(ga) ga[i] = g[i]; }
};

// v{*,a,b} x m{c,d,a,b} -> {*,c,d}
struct VecMatMulT : public DiffAB {
    VecMatMulT(PDiff a, PDiff b) : DiffAB(a, b) {}

    void Setup() override {
        auto as = a->shape(), bs = b->shape();
        Check(as.size() - 1 < bs.size());
        Check(as.pop_front() == bs.last(as.size() - 1));
        Reshape(bs.first(bs.size() - as.size() + 1));
    }

    void Forward() override {
        for (size_t i = 0; i < vb.size(); i += va.size()) {
            float s = 0;
            for (size_t j = 0; j < va.size(); j++) s += va[j] * vb[i + j];
            v[i] = s;
        }
    }

    void Backward() override {
        for (size_t j = 0; j < va.size(); j++) {
            float s = 0;
            for (size_t i = 0; i < vb.size(); i += va.size()) s += g[i] * vb[i + j];
            ga[j] += s;
        }
        //  TODO
    }
};

Declare2(VecMatMul, VecMatMulT);

struct SumT : public DiffA {
    SumT(PDiff a) : DiffA(a) {}
    void Setup() override { Reshape({1}); }

    void Forward() override {
        tensor::type sum = 0;
        EACH(va) sum += va[i];
        v[0] = sum;
    }

    void Backward() override { EACH(ga) ga[i] += g[0]; }
};

struct SizeT : public DiffA {
    SizeT(PDiff a) : DiffA(a) {}
    void Setup() override { v.reshape({1}); }
    void Forward() override { v[0] = va.size(); }
};

Declare1(Sum);
Declare1(Size);

#if 0
struct MeanT : public DiffA {
    MeanT(PDiff a) : DiffA(a) { }

    void Setup() override {
        Reshape({1});
    }

    void Forward() override {
        v[0] = ::Sum<tensor::type>(va) / va.size();
    }

    void Backward() override {
        EACH(ga) ga[i] += g[0] / va.size();
    }
};
#endif

// Declare1(Mean, MeanT)
// Declare1(SoftMax, SoftMaxT)

#if 0
template<bool mean>
struct SumSquareDiff : public DiffAB {
    SumSquareDiff(PDiff a, PDiff b) : DiffAB(a, b) { }

    void Setup() override {
        Check(a->shape() == b->shape() || b->size() == 1);
        Reshape({1});
    }

    void Forward() override {
        tensor::type sum = 0;
        if (vb.size() == 1)
            EACH(v) sum += sqr(va[i] - vb[0]);
        else
            EACH(v) sum += sqr(va[i] - vb[i]);
        if (mean) sum /= va.size();
        v[0] = sum;
    }

    void Backward() override {
        auto d = 2 * g[0];
        if (mean) d /= va.size();
        if (vb.size() == 1) {
            EACH(ga) ga[i] += (va[i] - vb[0]) * d;
            if (gb) for (auto i : range(va)) gb[0] += (vb[0] - va[i]) * d;
        } else {
            EACH(ga) ga[i] += (va[i] - vb[i]) * d;
            EACH(gb) gb[i] += (vb[i] - va[i]) * d;
        }
    }
};
#endif

#undef va
#undef vb
#undef vc

#undef ga
#undef gb
#undef gc

inline PDiff Mean(PDiff a) { return Sum(a) / Size(a); }

inline PDiff MeanSquareError(PDiff a, PDiff b) { return Sum(Sqr(a - b)) / Size(a); }

inline PDiff Softmax(PDiff a) {
    PDiff e = Exp(a);
    return e / Sum(e);
}

inline PDiff BatchNorm(PDiff a, double delta) {
    auto mean = Broadcast(Mean(a), a);
    auto stdev = Broadcast(Sqrt(delta + Mean(Sqr(a - mean))), a);
    return (a - mean) / stdev;
}

inline PDiff FullyConnected(PDiff a, uint16_t size, shared_ptr<Init> w_init = make_shared<Init>()) {
    auto w = Param(a->shape().remove_zeros().push_front(size), "w", w_init);
    auto b = Param({size}, "b");
    return VecMatMul(a, w) + b;
}

inline bool IsDiff(PDiff a) {
    const auto& e = *a.get();
    return typeid(e) == typeid(Diff);
}

inline bool IsData(PDiff a) { return IsDiff(a) && !a->g && a->shape().size() > 1 && a->shape()[0] == 0; }

inline bool IsParam(PDiff a) { return IsDiff(a) && a->g && a->shape().size() > 0 && a->shape()[0] != 0; }

inline bool IsConst(PDiff a) { return IsDiff(a) && !a->g && a->shape().size() > 0 && a->shape()[0] != 0; }

void Setup(span<PDiff> nodes);
void InitBatches(span<PDiff> nodes, int batch_size);
void ResetTicks(span<PDiff> nodes);
void Forward(span<PDiff> nodes);
void ResetGradients(span<PDiff> nodes);
void Backward(span<PDiff> nodes);
void GradientDescent(span<PDiff> nodes, const float alpha);
void CheckBounded(cspan<PDiff> nodes, const tensor::type limit);
void Print(cspan<PDiff> nodes);

using Metrics = unordered_map<string, float>;

struct Model {
    Model(PDiff loss, PDiff accuracy, int batch_size);

    PDiff loss;
    PDiff accuracy;
    vector<PDiff> nodes;
    // TODO filtered lists of nodes for forward, reset_gradients, backward, gradient_descent

    void Forward() { ::Forward(nodes); }
    void ResetGradients() { ::ResetGradients(nodes); }
    void Backward() { ::Backward(nodes); }
    void GradientDescent(tensor::type alpha) { ::GradientDescent(nodes, alpha); }
    void Print() const { ::Print(nodes); }

    Metrics Epoch(cspan<pair<PDiff, tensor>> data, const float alpha, const size_t seed);
};
