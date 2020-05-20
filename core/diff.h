#pragma once
#include <core/auto.h>
#include <core/range.h>
#include <core/std.h>
#include <core/tensor.h>

// Stages of Diff graph:
// - construction - final size of tensors is known
// - optimization
// - Forward() - v tensors populated
// - ResetGradients() - g tensors cleared
// - Backward() - g tensors accumulated
// - GradientDescent() - v tensors adjusted using g tensors

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
    virtual void Forward() {}
    virtual void Backward() {}

    void Reshape(tensor_shape s) {
        v.reshape(s);
        g.reshape(s);
    }

    string str() const { return format("%s %s", TypeName(*this), string(v.shape())); }
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

struct DiffSimpleA : public Diff {
    DiffSimpleA(PDiff a) : a(a) { }
    array<PDiff, 3> Inputs() override { return {a, nullptr, nullptr}; }

    PDiff a;
};

struct DiffA : public DiffSimpleA {
    DiffA(PDiff a) : DiffSimpleA(a) { Reshape(a->size()); }
};

struct DiffSimpleAB : public Diff {
    DiffSimpleAB(PDiff a, PDiff b) : a(a), b(b) { }
    array<PDiff, 3> Inputs() override { return {a, b, nullptr}; }

    PDiff a, b;
};

bool IsBroadcastable(tensor_shape a, tensor_shape b);
PDiff Broadcast(PDiff a, tensor_shape b);

struct DiffAB : public DiffSimpleAB {
    DiffAB(PDiff a, PDiff b) : DiffSimpleAB(a, b) {
        if (a->shape().remove_ones() != b->shape().remove_ones()) {
            println("%s %s %s %s\n", a->shape(), b->shape(), TypeName(*a), TypeName(*b));
            if (a->size() > b->size()) b = Broadcast(b, a->shape());
            if (a->size() < b->size()) a = Broadcast(a, b->shape());
            Check((a->size() == 1 && b->size() == 1) || a->shape() == b->shape(), format("%s %s %s", a->shape(), b->shape(), TypeName(*this)));
        }
        Reshape(a->shape());
    }
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
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) * g[i]; }
};

struct LeakyReluT : public DiffA {
    LeakyReluT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (va[i] / 10); }
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) ? g[i] : (g[i] / 10); }
};

struct LogisticT : public DiffA {
    LogisticT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = 1 / (1 + exp(-va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i] * (1 - v[i]); }
};

struct TanhT : public DiffA {
    TanhT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = tanh(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * (1 - sqr(v[i])); }
};

struct ErfT : public DiffA {
    ErfT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = erf(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * (2 / sqrt(PI)) * exp(-sqr(va[i])); }
};

struct AtanT : public DiffA {
    AtanT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = atan(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / (1 + sqr(va[i])); }
};

// x / (1 + abs(x))
struct SaxT : public DiffA {
    SaxT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = va[i] / (1 + abs(va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] / sqr(1 + abs(va[i])); }
};

template<typename T>
T cube(T a) { return a * a * a; }

// 1/sqrt(1+x^2)
struct SoxT : public DiffA {
    SoxT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = 1 / sqrt(1 + sqr(va[i])); }
    void Backward() override { EACH(ga) ga[i] -= g[i] * va[i] * cube(v[i]); }
};

// TODO non-linealities to try:
// atan(pi*x/2)*2/pi

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
    void Backward() override { EACH(ga) ga[i] += g[i] * sign(va[i]); }
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
Declare1(Logistic);
Declare1(Tanh);
Declare1(Erf);
Declare1(Atan);
Declare1(Sax);
Declare1(Sox);
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
    BroadcastFMA(PDiff a, PDiff b, PDiff c) : DiffABC(a, b, c) {
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

inline tensor::type Sum(const tensor a, tensor::type s = 0) {
    EACH(a) s += a[i];
    return s;
}

inline tensor::type Dot(const tensor a, const tensor b, tensor::type s = 0) {
    DCheck(a.size() == b.size(), "arguments");
    EACH(a) s += a[i] * b[i];
    return s;
}

struct AddVV : public DiffAB {
    AddVV(PDiff a, PDiff b) : DiffAB(a, b) { }
    void Forward() override { EACH(v) v[i] = va[i] + vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] += g[i];
    }
};

struct AddVS : public DiffSimpleAB {
    AddVS(PDiff a, PDiff b) : DiffSimpleAB(a, b) {
        Check(b->size() == 1);
        Reshape(a->shape());
    }
    void Forward() override { EACH(v) v[i] = va[i] + vb[0]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        if (gb) gb[0] = Sum(g, gb[0]);
    }
};

struct AddMV : public DiffSimpleAB {
    AddMV(PDiff a, PDiff b) : DiffSimpleAB(a, b) {
        Check(a->shape().pop_front() == b->shape());
        Reshape(a->shape());
    }
    void Forward() override {
        for (auto i : range<size_t>(0, va.size(), vb.size())) {
            for (auto j : range(vb.size())) v[i + j] = va[i + j] + vb[j];
        }
    }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        const auto m = va.shape()[0];
        EACH(gb) gb[i] += g[i] * m;
    }
};

inline PDiff operator+(PDiff a, PDiff b) {
    if (a->size() == 1) return make_shared<AddVS>(b, a);
    if (b->size() == 1) return make_shared<AddVS>(a, b);
    if (a->shape() == b->shape()) return make_shared<AddVV>(a, b);
    if (a->size() > b->size()) return make_shared<AddMV>(a, b);
    return make_shared<AddMV>(b, a);
}

inline auto operator+(tensor::type a, PDiff b) { return make_shared<AddVS>(b, Const(a)); }
inline auto operator+(PDiff a, tensor::type b) { return make_shared<AddVS>(a, Const(b)); }

struct Sub : public DiffAB {
    Sub(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] - vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] += -g[i];
    }
};

inline auto operator-(PDiff a, PDiff b) {
    return make_shared<Sub>(a, b);
}

struct MulVS : public DiffSimpleAB {
    MulVS(PDiff a, PDiff b) : DiffSimpleAB(a, b) {
        Check(b->size() == 1);
        Reshape(a->shape());
    }
    void Forward() override { EACH(v) v[i] = va[i] * vb[0]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i] * vb[0];
        if (gb) gb[0] = Dot(g, va, gb[0]);
    }
};

struct MulVV : public DiffAB {
    MulVV(PDiff a, PDiff b) : DiffAB(a, b) { }
    void Forward() override { EACH(v) v[i] = va[i] * vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i] * vb[i];
        EACH(gb) gb[i] += g[i] * va[i];
    }
};

inline PDiff operator*(PDiff a, PDiff b) {
    if (a->size() == 1) return make_shared<MulVS>(b, a);
    if (b->size() == 1) return make_shared<MulVS>(a, b);
    return make_shared<MulVV>(a, b);
}

inline auto operator*(tensor::type a, PDiff b) { return make_shared<MulVS>(b, Const(a)); }
inline auto operator*(PDiff a, tensor::type b) { return make_shared<MulVS>(a, Const(b)); }

struct Div : public DiffAB {
    Div(PDiff a, PDiff b) : DiffAB(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] / vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i] / vb[i];
        EACH(gb) gb[i] -= g[i] / (va[i] * va[i]);
    }
};

inline auto operator/(PDiff a, PDiff b) {
    return make_shared<Div>(a, b);
}

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

#define CONST_OVERLOAD(OP)                                                     \
    inline auto operator OP(tensor::type a, PDiff b) { return Broadcast(Const(a), b->shape()) OP b; } \
    inline auto operator OP(PDiff a, tensor::type b) { return a OP Broadcast(Const(b), a->shape()); }

CONST_OVERLOAD(-);
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

struct Concat : public DiffSimpleAB {
    Concat(PDiff a, PDiff b) : DiffSimpleAB(a, b) {
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

struct Conv1D : public DiffSimpleAB {
    Conv1D(PDiff a, PDiff b, int offset) : DiffSimpleAB(a, b), offset(offset) {
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

struct Deconv1D : public DiffSimpleAB {
    Deconv1D(PDiff a, PDiff b) : DiffSimpleAB(a, b) {}
};

struct MaxPool1D : public DiffSimpleA {
    MaxPool1D(PDiff a) : DiffSimpleA(a) {
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

struct MaxPool2D : public DiffSimpleA {
    MaxPool2D(PDiff a);
    void Forward() override;
    void Backward() override;
};

struct Reshape : public DiffSimpleA {
    Reshape(PDiff a, tensor_shape shape) : DiffSimpleA(a) {
        v.reshape(shape);
        Check(a->size() == v.size());
    }

    void Forward() override { EACH(v) v[i] = va[i]; }
    void Backward() override { EACH(ga) ga[i] = g[i]; }
};

// m can be zero, one or more dimensions
// p and q are exactly one dimension
// v{m,p} x m{q,p} -> {m,q}
struct VecMatMulT : public DiffSimpleAB {
    VecMatMulT(PDiff a, PDiff b) : DiffSimpleAB(a, b) {
        auto as = a->shape(), bs = b->shape();
        Check(as.size() > 0);
        Check(bs.size() == 2);
        Check(as.back() == bs.back());
        Reshape(as.pop_back().push_back(bs[0]));
    }

    void Forward() override {
        for (size_t i = 0; i < vb.size(); i += va.size()) {
            tensor::type s = 0;
            for (size_t j = 0; j < va.size(); j++) s += va[j] * vb[i + j];
            v[i] = s;
        }
    }

    void Backward() override {
        EACH (ga) {
            tensor::type s = 0;
            for (size_t j = 0; j < vb.size(); j += va.size()) s += g[j] * vb[j + i];
            ga[i] += s;
        }
        if (gb) {
            //  TODO
        }
    }
};

Declare2(VecMatMul, VecMatMulT);

struct SumT : public DiffSimpleA {
    SumT(PDiff a) : DiffSimpleA(a) { Reshape({1}); }
    void Forward() override { v[0] = Sum(va); }
    void Backward() override { EACH(ga) ga[i] += g[0]; }
};

Declare1(Sum);

#if 0
struct MeanT : public DiffSimpleA {
    MeanT(PDiff a) : DiffSimpleA(a) { Reshape({1}); }
    void Forward() override { v[0] = Sum(va) / va.size(); }
    void Backward() override { EACH(ga) ga[i] += g[0] / va.size(); }
};
#endif

// Declare1(Mean, MeanT)
// Declare1(SoftMax, SoftMaxT)

#if 0
template<bool mean>
struct SumSquareDiff : public DiffSimpleAB {
    SumSquareDiff(PDiff a, PDiff b) : DiffSimpleAB(a, b) {
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

inline PDiff Mean(PDiff a) { return Sum(a) / a->size(); }

inline PDiff MeanSquareError(PDiff a, PDiff b) { return Sum(Sqr(a - b)) / a->size(); }

inline PDiff Softmax(PDiff a) {
    PDiff e = Exp(a);
    return e / Sum(e);
}

inline PDiff BatchNorm(PDiff a, double delta) {
    auto mean = Broadcast(Mean(a), a->shape());
    auto stdev = Broadcast(Sqrt(delta + Mean(Sqr(a - mean))), a->shape());
    return (a - mean) / stdev;
}

inline PDiff FullyConnected(PDiff a, uint16_t size, shared_ptr<Init> w_init = make_shared<Init>()) {
    auto w = Param(tensor_shape(size, a->shape().back()), "w", w_init);
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
    Model(PDiff loss, PDiff accuracy);

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
