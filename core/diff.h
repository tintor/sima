#pragma once
#include <core/auto.h>
#include <core/range.h>
#include <core/std.h>
#include <core/tensor.h>

// Stages of Diff graph:
// - construction - final size of tensors is known
// - optimization
// - Forward(bool) - v tensors populated
// - ResetGradients() - g tensors cleared
// - Backward() - g tensors accumulated
// - Optimize() - v tensors adjusted using g tensors

struct Diff;

using PDiff = std::shared_ptr<Diff>;

struct Diff {
    virtual vector<PDiff> Inputs() { return {}; }

    static thread_local bool has_overload;
    virtual void BeginEpoch(bool training) { has_overload = false; }
    virtual void Forward(bool training) { has_overload = false; }
    virtual void Backward() { has_overload = false; }
    virtual void EndEpoch(bool training) { has_overload = false; }

    void SetBatchSize(uint batch) {
        if (batched()) {
            const dim4 s = v.shape().set(0, batch, 'b');
            v.reshape(s);
            if (g) g.reshape(s);
        }
    }

    void Reshape(dim4 s) {
        v.reshape(s);
        g.reshape(s);
    }

    operator string() const { return format("%s %s", TypeName(*this), v.shape().str()); }

    dim_t dim(uint i) const { return v.shape()[i]; }
    auto shape() const { return v.shape(); }
    auto ndims() const { return v.ndims(); }
    auto elements() const { return v.elements(); }
    bool batched() const { return ndims() > 0 && shape().name(0) == 'b'; }

    bool inference = true;
    string name;
    vtensor v, g;

    ulong forward_ticks = 0;
    ulong backward_ticks = 0;
};

inline PDiff operator<<(PDiff a, string_view name) {
    a->name = name;
    return a;
}

#define FOR(i, END) for (auto i : range(END))
#define EACH(V) FOR(i, V.elements())

#define Declare1(Func) \
    inline PDiff Func(PDiff a) { return make_shared<Func##T>(a); }

struct Diff1 : public Diff {
    Diff1(PDiff a) : a(a), va(a->v), ga(a->g) {}
    vector<PDiff> Inputs() override { return {a}; }
    PDiff a;
    vtensor& va;
    vtensor& ga;
};

struct Diff2 : public Diff1 {
    Diff2(PDiff a, PDiff b) : Diff1(a), b(b), vb(b->v), gb(b->g) {}
    vector<PDiff> Inputs() override { return {a, b}; }
    PDiff b;
    vtensor& vb;
    vtensor& gb;
};

struct Diff3 : public Diff2 {
    Diff3(PDiff a, PDiff b, PDiff c) : Diff2(a, b), c(c), vc(c->v), gc(c->g) {}
    vector<PDiff> Inputs() override { return {a, b, c}; }
    PDiff c;
    vtensor& vc;
    vtensor& gc;
};

struct DiffA : public Diff1 {
    DiffA(PDiff a) : Diff1(a) { Reshape(a->shape()); }
};

struct Diff_vv : public Diff2 {
    Diff_vv(PDiff a, PDiff b) : Diff2(a, b) {
        Check(a->shape() == b->shape(), format("%s vs %s", string(a->shape()), string(b->shape())));
        Reshape(a->shape());
    }
};

struct Diff_sv : public Diff2 {
    Diff_sv(PDiff a, PDiff b) : Diff2(a, b) {
        Check(a->elements() == 1, string(a->shape()));
        Reshape(b->shape());
    }
};

struct Diff_vs : public Diff2 {
    Diff_vs(PDiff a, PDiff b) : Diff2(a, b) {
        Check(b->elements() == 1, string(b->shape()));
        Reshape(a->shape());
    }
};

bool IsBroadcastable(dim4 a, dim4 b);
PDiff Broadcast(PDiff a, dim4 b);

// Ignores gradients and batches.
inline PDiff Const(tensor::type c) {
    auto p = make_shared<Diff>();
    p->v.reshape({});
    p->v[0] = c;
    return p;
}

inline PDiff Const(dim4 shape, tensor::type c) {
    auto p = make_shared<Diff>();
    p->v.reshape(shape, c);
    return p;
}

inline PDiff Const(const tensor c) {
    auto p = make_shared<Diff>();
    p->v = c;
    return p;
}

inline bool IsConst(PDiff a) { return IsType<Diff>(a) && !a->g; }

struct Init {
    virtual tensor::type get() { return 0; }
};

struct NormalInit : public Init {
    NormalInit(tensor::type variance, std::mt19937_64& random) : dis(0, variance), random(random) {}
    tensor::type get() override { return dis(random); }

    std::normal_distribution<tensor::type> dis;
    std::mt19937_64& random;
};

struct UniformInit : public Init {
    UniformInit(tensor::type min, tensor::type max, std::mt19937_64& random) : dis(min, max), random(random) {}
    tensor::type get() override { return dis(random); }

    std::uniform_real_distribution<tensor::type> dis;
    std::mt19937_64& random;
};

struct ParamT : public Diff {};

// Accepts gradients, ignores batches. Learnable parameter, can be saved.
inline PDiff Param(dim4 shape, shared_ptr<Init> init = nullptr) {
    auto p = make_shared<ParamT>();
    p->Reshape(shape);
    if (init) EACH(p->v) p->v[i] = init->get();
    return p;
}

inline PDiff Param(dim4 shape, tensor::type init) {
    auto p = make_shared<ParamT>();
    p->v.reshape(shape, init);
    p->g.reshape(shape);
    return p;
}

// Input and reference value. Not saved.
// Ignores gradients, requires batches.
inline PDiff Data(dim4 shape) {
    auto p = make_shared<Diff>();
    p->v.reshape(shape.push_front(1, 'b'));
    return p;
}

struct GaussianT : public Diff {
    GaussianT(dim4 shape, tensor::type mean, tensor::type stdev, size_t seed) : normal(mean, stdev), random(seed) {
        v.reshape(shape);
    }
    void Forward(bool) override { EACH(v) v[i] = normal(random); }
    std::normal_distribution<tensor::type> normal;
    std::mt19937_64 random;
};

inline PDiff Gaussian(dim4 shape, tensor::type mean, tensor::type stdev, size_t seed) {
    return make_shared<GaussianT>(shape, mean, stdev, seed);
}

struct ReluT : public DiffA {
    ReluT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = (va[i] > 0) ? va[i] : 0; }
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) * g[i]; }
};
Declare1(Relu);

PDiff operator+(PDiff a, PDiff b);
inline PDiff NoisyRelu(PDiff a, size_t seed) { return Relu(a + Gaussian(a->shape(), 0, 1, seed)); }

// SmoothRelu
struct SoftplusT : public DiffA {
    SoftplusT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = log(1 + std::exp(va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] / (1 + exp(-va[i])); }
};
Declare1(Softplus);

struct LeakyReluT : public DiffA {
    LeakyReluT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (va[i] / 100); }
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) ? g[i] : (g[i] / 100); }
};
Declare1(LeakyRelu);

struct ParametricReluT : public Diff_vv {
    ParametricReluT(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward(bool) override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (va[i] * vb[i]); }
    void Backward() override {
        EACH(ga) ga[i] += (va[i] > 0) ? g[i] : (g[i] * vb[i]);
        EACH(gb) gb[i] += (va[i] > 0) ? 0 : (g[i] * va[i]);
    }
};
inline PDiff ParametricRelu(PDiff a, PDiff b) { return make_shared<ParametricReluT>(a, b); }

struct ELUT : public DiffA {
    ELUT(PDiff a, double k) : DiffA(a), k(k) {}
    void Forward(bool) override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (k * (std::exp(va[i] - 1))); }
    void Backward() override { Check(false, "not implemented"); }
    double k;
};
inline PDiff ELU(PDiff a, double k) { return make_shared<ELUT>(a, k); }

struct LogisticT : public DiffA {
    LogisticT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = 1 / (1 + std::exp(-va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i] * (1 - v[i]); }
};
Declare1(Logistic);

struct ClampLogisticT : public DiffA {
    ClampLogisticT(PDiff a, tensor::type limit) : DiffA(a), limit(limit) {}
    void Forward(bool) override { EACH(v) v[i] = 1 / (1 + std::exp(std::clamp(-va[i], -limit, limit))); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i] * (1 - v[i]); }
    tensor::type limit;
};

inline PDiff Logistic(PDiff a, float limit) { return make_shared<ClampLogisticT>(a, limit); }

struct TanhT : public DiffA {
    TanhT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::tanh(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * (1 - sqr(v[i])); }
};
Declare1(Tanh);

struct ErfT : public DiffA {
    ErfT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::erf(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * tensor::type(2 / sqrt(PI)) * std::exp(-sqr(va[i])); }
};
Declare1(Erf);

// TODO variation: atan(pi*x/2)*2/pi
struct AtanT : public DiffA {
    AtanT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::atan(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / (1 + sqr(va[i])); }
};
Declare1(Atan);

// x / (1 + abs(x))
struct SaxT : public DiffA {
    SaxT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = va[i] / (1 + std::abs(va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] / sqr(1 + std::abs(va[i])); }
};
Declare1(Sax);

// 1/sqrt(1+x^2)
struct SoxT : public DiffA {
    SoxT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = 1 / std::sqrt(1 + sqr(va[i])); }
    void Backward() override { EACH(ga) ga[i] -= g[i] * va[i] * cube(v[i]); }
};
Declare1(Sox);

struct SqrT : public DiffA {
    SqrT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = va[i] * va[i]; }
    void Backward() override { EACH(ga) ga[i] += g[i] * 2 * va[i]; }
};
Declare1(Sqr);

struct SqrtT : public DiffA {
    SqrtT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::sqrt(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / 2 / v[i]; }
};
Declare1(Sqrt);

struct ExpT : public DiffA {
    ExpT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::exp(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i]; }
};
Declare1(Exp);

struct LogT : public DiffA {
    LogT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::log(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / va[i]; }
};
Declare1(Log);

PDiff Pow(PDiff a, PDiff b);

struct AbsT : public DiffA {
    AbsT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::abs(va[i]); }
    void Backward() override { EACH(ga) ga[i] += std::copysign(g[i], va[i]); }
};
Declare1(Abs);

struct CosT : public DiffA {
    CosT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::cos(va[i]); }
    void Backward() override { EACH(ga) ga[i] -= std::sin(va[i]); }
};
Declare1(Cos);

struct SinT : public DiffA {
    SinT(PDiff a) : DiffA(a) {}
    void Forward(bool) override { EACH(v) v[i] = std::sin(va[i]); }
    void Backward() override { EACH(ga) ga[i] += std::cos(va[i]); }
};
Declare1(Sin);

inline tensor::type Sum(const tensor a, tensor::type s = 0) {
    EACH(a) s += a[i];
    return s;
}

inline tensor::type Dot(const tensor a, const tensor b, tensor::type s = 0) {
    DCheck(a.size() == b.size(), "arguments");
    EACH(a) s += a[i] * b[i];
    return s;
}

PDiff operator+(PDiff a, PDiff b);
PDiff operator+(tensor::type a, PDiff b);
PDiff operator+(PDiff a, tensor::type b);

PDiff operator-(PDiff a, PDiff b);
PDiff operator-(PDiff a);
PDiff operator-(tensor::type a, PDiff b);
PDiff operator-(PDiff a, tensor::type b);

PDiff operator*(PDiff a, PDiff b);
PDiff operator*(tensor::type a, PDiff b);
PDiff operator*(PDiff a, tensor::type b);

PDiff Inv(PDiff a);
PDiff operator/(PDiff a, PDiff b);
PDiff operator/(tensor::type a, PDiff b);
PDiff operator/(PDiff a, tensor::type b);

inline PDiff& operator+=(PDiff& a, PDiff b) { return a = a + b; }
inline PDiff& operator-=(PDiff& a, PDiff b) { return a = a - b; }
inline PDiff& operator*=(PDiff& a, PDiff b) { return a = a * b; }
inline PDiff& operator/=(PDiff& a, PDiff b) { return a = a / b; }

inline PDiff& operator+=(PDiff& a, tensor::type b) { return a = a + b; }
inline PDiff& operator-=(PDiff& a, tensor::type b) { return a = a - b; }
inline PDiff& operator*=(PDiff& a, tensor::type b) { return a = a * b; }
inline PDiff& operator/=(PDiff& a, tensor::type b) { return a = a / b; }

#define RELATION(NAME, OP)                                                                            \
    struct NAME : public Diff_vv {                                                                    \
        NAME(PDiff a, PDiff b) : Diff_vv(a, b) {}                                                     \
        void Forward(bool) override { EACH(v) v[i] = (va[i] OP vb[i]) ? 1 : 0; }                      \
    };                                                                                                \
    inline PDiff operator OP(PDiff a, PDiff b) { return make_shared<NAME>(a, b); }                    \
    inline auto operator OP(tensor::type a, PDiff b) { return Broadcast(Const(a), b->shape()) OP b; } \
    inline auto operator OP(PDiff a, tensor::type b) { return a OP Broadcast(Const(b), a->shape()); }

RELATION(Greater, >);
RELATION(Less_, <);
RELATION(GreaterOrEqual, >=);
RELATION(LessOrEqual, <=);

#undef RELATION

PDiff Min(PDiff a, PDiff b);
inline auto Min(tensor::type a, PDiff b) { return Min(Broadcast(Const(a), b->shape()), b); }
inline auto Min(PDiff a, tensor::type b) { return Min(a, Broadcast(Const(b), a->shape())); }

PDiff Max(PDiff a, PDiff b);
inline auto Max(tensor::type a, PDiff b) { return Max(Broadcast(Const(a), b->shape()), b); }
inline auto Max(PDiff a, tensor::type b) { return Max(a, Broadcast(Const(b), a->shape())); }

PDiff Concat(PDiff a, PDiff b);
PDiff Reshape(PDiff a, dim4 shape);
PDiff VecMatMul(PDiff a, PDiff b);

PDiff Sum(PDiff a);
PDiff Mean(PDiff a);
PDiff Stdev(PDiff a, PDiff mean_a, tensor::type k = 0);
inline PDiff Stdev(PDiff a) { return Stdev(a, Mean(a)); }

// TODO Unary horizontal Min and Max!

PDiff RunningAverage(PDiff a, float k);
PDiff EpochMean(PDiff a, float init);
PDiff MeanSquareError(PDiff a, PDiff b);

PDiff ValueCmp(PDiff a, PDiff b);
pair<PDiff, PDiff> GradCmp(PDiff a);

PDiff BinaryCrossEntropy(PDiff ref, PDiff out);
PDiff BinaryAccuracy(PDiff ref, PDiff out);

inline PDiff Softmax(PDiff a) {
    PDiff e = Exp(a);
    return e / Sum(e) << "softmax";
}

PDiff BatchNorm(PDiff a, tensor::type k);

inline PDiff Bias(PDiff a) {
    const dim4 s = a->batched() ? a->shape().pop_front() : a->shape();
    return a + Param(s) << "bias";
}

inline PDiff Affine(PDiff a, uint channels, shared_ptr<Init> init) {
    auto w = Param({channels, a->shape().back()}, init);
    return VecMatMul(a, w) << "affine";
}

inline PDiff FullyConnected(PDiff a, uint channels, shared_ptr<Init> init) {
    return Bias(Affine(a, channels, init));
}

float ComputePolynomial(float x, cspan<float> poly);
float ComputePolynomialDeriv(float x, cspan<float> poly);

PDiff Kronecker(PDiff a);
PDiff Polynomial(PDiff a, uint degree, uint channels, shared_ptr<Init> init);
