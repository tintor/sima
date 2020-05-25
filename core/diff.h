#pragma once
#include <core/auto.h>
#include <core/property.h>
#include <core/range.h>
#include <core/std.h>
#include <core/tensor.h>
#include <arrayfire.h>

// Stages of Diff graph:
// - construction - final size of tensors is known
// - optimization
// - Forward() - v tensors populated
// - ResetGradients() - g tensors cleared
// - Backward() - g tensors accumulated
// - Optimize() - v tensors adjusted using g tensors

struct Diff;

using PDiff = std::shared_ptr<Diff>;

struct Diff;

struct Diff {
    virtual vector<PDiff> Inputs() { return {}; }

    static thread_local bool has_overload;
    virtual void Forward() { has_overload = false; }
    virtual void Backward() { has_overload = false; }
    virtual void EndEpoch() { has_overload = false; }

    void Reshape(tensor_shape s) {
        v.reshape(s);
        g.reshape(s);
    }

    operator string() const { return format("%s %s", TypeName(*this), string(v.shape())); }

    TProperty(Shape, Diff) {
        operator auto() const { return parent->v.shape(); }
        auto operator->() const { return &parent->v.shape(); }
        auto operator()(int i) const { return parent->v.shape()[i]; }
        bool operator==(const Shape& o) const { return parent->v.shape() == o; }
        bool operator==(const tensor_shape& o) const { return parent->v.shape() == o; }
        operator string() const { return string(parent->v.shape()); }
    }
    shape;

    Property(Diff) {
        operator uint() const { return parent->v.shape().size; }
    }
    rank;

    Property(Diff) {
        operator size_t() const { return parent->v.size; }
    }
    size;

    string name;
    vtensor v, g;

    ulong forward_ticks = 0;
    ulong backward_ticks = 0;
};

inline PDiff operator<<(PDiff a, string_view name) {
    a->name = name;
    return a;
}

#define EACH(V) for (auto i : range(size_t(V.size)))

#define Declare1(Func) \
    inline PDiff Func(PDiff a) { return make_shared<Func##T>(a); }
#define Declare2(Func, Type) \
    inline PDiff Func(PDiff a, PDiff b) { return make_shared<Type>(a, b); }

#define TensorProperty(TENSOR, TYPE, PARENT)                                               \
    TProperty(TYPE, PARENT) {                                                              \
        tensor::type operator[](size_t i) const { return parent->TENSOR[i]; }              \
        tensor::type& operator[](size_t i) { return parent->TENSOR[i]; }                   \
        tensor::type operator()(size_t i, size_t j) const { return parent->TENSOR(i, j); } \
        tensor::type& operator()(size_t i, size_t j) { return parent->TENSOR(i, j); }      \
        Property(TYPE) {                                                                   \
            operator size_t() const { return parent->parent->TENSOR.size; }                \
        }                                                                                  \
        size;                                                                              \
        const auto& shape() const { return parent->TENSOR.shape(); }                       \
        operator bool() const { return parent->TENSOR; }                                   \
        operator tensor() { return parent->TENSOR; }                                       \
        operator vtensor&() { return parent->TENSOR; }                                     \
    }

struct Diff1 : public Diff {
    Diff1(PDiff a) : a(a) {}
    vector<PDiff> Inputs() override { return {a}; }
    PDiff a;

    TensorProperty(a->v, AV, Diff1) va;
    TensorProperty(a->g, AG, Diff1) ga;
};

struct Diff2 : public Diff1 {
    Diff2(PDiff a, PDiff b) : Diff1(a), b(b) {}
    vector<PDiff> Inputs() override { return {a, b}; }
    PDiff b;

    TensorProperty(b->v, BV, Diff2) vb;
    TensorProperty(b->g, BG, Diff2) gb;
};

struct Diff3 : public Diff2 {
    Diff3(PDiff a, PDiff b, PDiff c) : Diff2(a, b), c(c) {}
    vector<PDiff> Inputs() override { return {a, b, c}; }
    PDiff c;

    TensorProperty(b->v, CV, Diff2) vc;
    TensorProperty(b->g, CG, Diff2) gc;
};

struct DiffA : public Diff1 {
    DiffA(PDiff a) : Diff1(a) { Reshape(a->shape); }
};

struct Diff_vv : public Diff2 {
    Diff_vv(PDiff a, PDiff b) : Diff2(a, b) {
        Check(a->shape == b->shape);
        Reshape(a->shape);
    }
};

struct Diff_sv : public Diff2 {
    Diff_sv(PDiff a, PDiff b) : Diff2(a, b) {
        Check(a->size == 1);
        Reshape(b->shape);
    }
};

struct Diff_vs : public Diff2 {
    Diff_vs(PDiff a, PDiff b) : Diff2(a, b) {
        Check(b->size == 1);
        Reshape(a->shape);
    }
};

bool IsBroadcastable(tensor_shape a, tensor_shape b);
PDiff Broadcast(PDiff a, tensor_shape b);

// Ignores gradients and batches.
inline PDiff Const(tensor::type c) {
    auto p = make_shared<Diff>();
    p->v.reshape({1});
    p->v[0] = c;
    return p;
}

inline PDiff Const(tensor_shape shape, tensor::type c) {
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

struct ParamT : public Diff {
};

// Accepts gradients, ignores batches. Learnable parameter, can be saved.
inline PDiff Param(tensor_shape shape, shared_ptr<Init> init = nullptr) {
    auto p = make_shared<ParamT>();
    Check(shape.size > 0, "Parameter shape.size() must be non-zero");
    p->Reshape(shape);
    if (init) EACH(p->v) p->v[i] = init->get();
    return p;
}

inline PDiff Param(tensor_shape shape, tensor::type init) {
    auto p = make_shared<ParamT>();
    Check(shape.size > 0, "Parameter shape.size() must be non-zero");
    p->v.reshape(shape, init);
    p->g.reshape(shape);
    return p;
}

// Input and reference value. Not saved.
// Ignores gradients, requires batches.
inline PDiff Data(tensor_shape shape) {
    auto p = make_shared<Diff>();
    p->v.reshape(shape);
    return p;
}

struct GaussianT : public Diff {
    GaussianT(tensor_shape shape, tensor::type mean, tensor::type stdev, size_t seed)
        : normal(mean, stdev), random(seed) {
        v.reshape(shape);
    }
    void Forward() override { EACH(v) v[i] = normal(random); }
    std::normal_distribution<tensor::type> normal;
    std::mt19937_64 random;
};

inline PDiff Gaussian(tensor_shape shape, tensor::type mean, tensor::type stdev, size_t seed) {
    return make_shared<GaussianT>(shape, mean, stdev, seed);
}

struct ReluT : public DiffA {
    ReluT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = (va[i] > 0) ? va[i] : 0; }
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) * g[i]; }
};
Declare1(Relu);

PDiff operator+(PDiff a, PDiff b);
inline PDiff NoisyRelu(PDiff a, size_t seed) { return Relu(a + Gaussian(a->shape, 0, 1, seed)); }

// SmoothRelu
struct SoftplusT : public DiffA {
    SoftplusT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = log(1 + std::exp(va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] / (1 + exp(-va[i])); }
};
Declare1(Softplus);

struct LeakyReluT : public DiffA {
    LeakyReluT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (va[i] / 100); }
    void Backward() override { EACH(ga) ga[i] += (va[i] > 0) ? g[i] : (g[i] / 100); }
};
Declare1(LeakyRelu);

struct ParametricReluT : public Diff_vv {
    ParametricReluT(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (va[i] * vb[i]); }
    void Backward() override {
        EACH(ga) ga[i] += (va[i] > 0) ? g[i] : (g[i] * vb[i]);
        EACH(gb) gb[i] += (va[i] > 0) ? 0 : (g[i] * va[i]);
    }
};
Declare2(ParametricRelu, ParametricReluT);

struct ELUT : public DiffA {
    ELUT(PDiff a, double k) : DiffA(a), k(k) {}
    void Forward() override { EACH(v) v[i] = (va[i] > 0) ? va[i] : (k * (std::exp(va[i] - 1))); }
    void Backward() override { Check(false, "not implemented"); }
    double k;
};
inline PDiff ELU(PDiff a, double k) { return make_shared<ELUT>(a, k); }

struct LogisticT : public DiffA {
    LogisticT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = 1 / (1 + std::exp(-va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i] * (1 - v[i]); }
};
Declare1(Logistic);

struct ClampLogisticT : public DiffA {
    ClampLogisticT(PDiff a, tensor::type limit) : DiffA(a), limit(limit) {}
    void Forward() override { EACH(v) v[i] = 1 / (1 + std::exp(std::clamp(-va[i], -limit, limit))); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i] * (1 - v[i]); }
    tensor::type limit;
};

inline PDiff Logistic(PDiff a, float limit) { return make_shared<ClampLogisticT>(a, limit); }

struct TanhT : public DiffA {
    TanhT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::tanh(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * (1 - sqr(v[i])); }
};
Declare1(Tanh);

struct ErfT : public DiffA {
    ErfT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::erf(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * tensor::type(2 / sqrt(PI)) * std::exp(-sqr(va[i])); }
};
Declare1(Erf);

// TODO variation: atan(pi*x/2)*2/pi
struct AtanT : public DiffA {
    AtanT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::atan(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / (1 + sqr(va[i])); }
};
Declare1(Atan);

// x / (1 + abs(x))
struct SaxT : public DiffA {
    SaxT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = va[i] / (1 + std::abs(va[i])); }
    void Backward() override { EACH(ga) ga[i] += g[i] / sqr(1 + std::abs(va[i])); }
};
Declare1(Sax);

// 1/sqrt(1+x^2)
struct SoxT : public DiffA {
    SoxT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = 1 / std::sqrt(1 + sqr(va[i])); }
    void Backward() override { EACH(ga) ga[i] -= g[i] * va[i] * cube(v[i]); }
};
Declare1(Sox);

struct SqrT : public DiffA {
    SqrT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = va[i] * va[i]; }
    void Backward() override { EACH(ga) ga[i] += g[i] * 2 * va[i]; }
};
Declare1(Sqr);

struct SqrtT : public DiffA {
    SqrtT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::sqrt(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / 2 / v[i]; }
};
Declare1(Sqrt);

struct ExpT : public DiffA {
    ExpT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::exp(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] * v[i]; }
};
Declare1(Exp);

struct LogT : public DiffA {
    LogT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::log(va[i]); }
    void Backward() override { EACH(ga) ga[i] += g[i] / va[i]; }
};
Declare1(Log);

struct AbsT : public DiffA {
    AbsT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::abs(va[i]); }
    void Backward() override { EACH(ga) ga[i] += std::copysign(g[i], va[i]); }
};
Declare1(Abs);

struct CosT : public DiffA {
    CosT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::cos(va[i]); }
    void Backward() override { EACH(ga) ga[i] -= std::sin(va[i]); }
};
Declare1(Cos);

struct SinT : public DiffA {
    SinT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = std::sin(va[i]); }
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

struct Add_vv : public Diff_vv {
    Add_vv(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] + vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] += g[i];
    }
};

struct Add_vs : public Diff_vs {
    Add_vs(PDiff a, PDiff b) : Diff_vs(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] + vb[0]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        if (gb) gb[0] += Sum(g, gb[0]);
    }
};

struct Add_mv : public Diff2 {
    Add_mv(PDiff a, PDiff b) : Diff2(a, b) {
        Check(a->shape->pop_front() == b->shape);
        Reshape(a->shape);
    }
    void Forward() override {
        for (auto i : range<size_t>(0, va.size, vb.size)) {
            for (auto j : range<size_t>(vb.size)) v[i + j] = va[i + j] + vb[j];
        }
    }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        const auto m = va.shape()[0];
        EACH(gb) gb[i] += g[i] * m;
    }
};

inline PDiff operator+(PDiff a, PDiff b) {
    if (a->size == 1 && b->size == 1) return make_shared<Add_vv>(a, b);
    if (a->size == 1) return make_shared<Add_vs>(b, a);
    if (b->size == 1) return make_shared<Add_vs>(a, b);
    if (a->shape == b->shape) return make_shared<Add_vv>(a, b);
    if (a->size > b->size) return make_shared<Add_mv>(a, b);
    return make_shared<Add_mv>(b, a);
}

inline auto operator+(tensor::type a, PDiff b) { return (a == 0) ? b : (b + Const(a)); }
inline auto operator+(PDiff a, tensor::type b) { return b + a; }

struct Sub_vv : public Diff_vv {
    Sub_vv(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] - vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] -= g[i];
    }
};

struct Sub_vs : public Diff_vs {
    Sub_vs(PDiff a, PDiff b) : Diff_vs(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] - vb[0]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i];
        if (gb) gb[0] -= Sum(g, gb[0]);
    }
};

struct Sub_sv : public Diff_sv {
    Sub_sv(PDiff a, PDiff b) : Diff_sv(a, b) {}
    void Forward() override { EACH(v) v[i] = va[0] - vb[i]; }
    void Backward() override {
        if (ga) ga[0] += Sum(g, ga[0]);
        EACH(gb) gb[i] -= g[i];
    }
};

struct Neg : public DiffA {
    Neg(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = -va[i]; }
    void Backward() override { EACH(ga) ga[i] -= g[i]; }
};

inline PDiff operator-(PDiff a, PDiff b) {
    if (a->size == 1 && b->size == 1) return make_shared<Sub_vv>(a, b);
    if (a->size == 1) return make_shared<Sub_sv>(a, b);
    if (b->size == 1) return make_shared<Sub_vs>(a, b);
    return make_shared<Sub_vv>(a, b);
}

inline PDiff operator-(tensor::type a, PDiff b) { return (a == 0) ? make_shared<Neg>(b) : (Const(a) - b); }
inline PDiff operator-(PDiff a, tensor::type b) { return (b == 0) ? a : (a - Const(b)); }
inline PDiff operator-(PDiff a) { return make_shared<Neg>(a); }

struct Mul_vs : public Diff_vs {
    Mul_vs(PDiff a, PDiff b) : Diff_vs(a, b) {
        Check(b->size == 1);
        Reshape(a->shape);
    }
    void Forward() override { EACH(v) v[i] = va[i] * vb[0]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i] * vb[0];
        if (gb) gb[0] += Dot(g, va);
    }
};

struct Mul_vv : public Diff_vv {
    Mul_vv(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] * vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i] * vb[i];
        EACH(gb) gb[i] += g[i] * va[i];
    }
};

inline PDiff operator*(PDiff a, PDiff b) {
    if (a->size == 1 && b->size == 1) return make_shared<Mul_vv>(a, b);
    if (a->size == 1) return make_shared<Mul_vs>(b, a);
    if (b->size == 1) return make_shared<Mul_vs>(a, b);
    return make_shared<Mul_vv>(a, b);
}

inline PDiff operator*(tensor::type a, PDiff b) {
    if (a == 1) return b;
    if (a == 0) return Const(b->shape, 0);
    return b * Const(a);
}

inline PDiff operator*(PDiff a, tensor::type b) {
    if (b == 1) return a;
    if (b == 0) return Const(a->shape, 0);
    return a * Const(b);
}

struct Div_vv : public Diff_vv {
    Div_vv(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] / vb[i]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i] / vb[i];
        EACH(gb) gb[i] -= g[i] * v[i] / vb[i];
    }
};

struct Div_vs : public Diff_vs {
    Div_vs(PDiff a, PDiff b) : Diff_vs(a, b) {}
    void Forward() override { EACH(v) v[i] = va[i] / vb[0]; }
    void Backward() override {
        EACH(ga) ga[i] += g[i] / vb[0];
        if (gb) EACH(g) gb[0] -= g[i] * v[i] / vb[0];
    }
};

struct Div_sv : public Diff_sv {
    Div_sv(PDiff a, PDiff b) : Diff_sv(a, b) {}
    void Forward() override { EACH(v) v[i] = va[0] / vb[i]; }
    void Backward() override {
        if (ga) EACH(g) ga[0] += g[i] / vb[i];
        EACH(gb) gb[i] -= g[i] * v[0] / sqr(vb[i]);
    }
};

struct InvT : public DiffA {
    InvT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = 1 / va[i]; }
    void Backward() override { EACH(ga) ga[i] -= g[i] * sqr(v[i]); }
};

Declare1(Inv);

inline PDiff operator/(PDiff a, PDiff b) {
    if (a->size == 1 && b->size == 1) return make_shared<Div_vv>(a, b);
    if (a->size == 1) return make_shared<Div_sv>(a, b);
    if (b->size == 1) return make_shared<Mul_vs>(a, Inv(b));
    return make_shared<Div_vv>(a, b);
}

inline auto operator/(tensor::type a, PDiff b) { return (a == 1) ? Inv(b) : (Const(a) / b); }
inline auto operator/(PDiff a, tensor::type b) { return a * Const(1 / b); }

inline PDiff& operator+=(PDiff& a, PDiff b) { return a = a + b; }
inline PDiff& operator-=(PDiff& a, PDiff b) { return a = a - b; }
inline PDiff& operator*=(PDiff& a, PDiff b) { return a = a * b; }
inline PDiff& operator/=(PDiff& a, PDiff b) { return a = a / b; }

inline PDiff& operator+=(PDiff& a, tensor::type b) { return a = a + b; }
inline PDiff& operator-=(PDiff& a, tensor::type b) { return a = a - b; }
inline PDiff& operator*=(PDiff& a, tensor::type b) { return a = a * b; }
inline PDiff& operator/=(PDiff& a, tensor::type b) { return a = a / b; }

#define RELATION(NAME, OP)                                                                          \
    struct NAME : public Diff_vv {                                                                  \
        NAME(PDiff a, PDiff b) : Diff_vv(a, b) {}                                                   \
        void Forward() override { EACH(v) v[i] = (va[i] OP vb[i]) ? 1 : 0; }                        \
    };                                                                                              \
    Declare2(operator OP, NAME);                                                                    \
    inline auto operator OP(tensor::type a, PDiff b) { return Broadcast(Const(a), b->shape) OP b; } \
    inline auto operator OP(PDiff a, tensor::type b) { return a OP Broadcast(Const(b), a->shape); }

RELATION(Greater, >);
RELATION(Less_, <);
RELATION(GreaterOrEqual, >=);
RELATION(LessOrEqual, <=);

#undef RELATION

struct MinT : public Diff_vv {
    MinT(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override { EACH(v) v[i] = min(va[i], vb[i]); }
    void Backward() override {
        EACH(ga) ga[i] += (va[i] < vb[i]) * g[i];
        EACH(gb) gb[i] += (vb[i] < va[i]) * g[i];
    }
};
Declare2(Min, MinT);
inline auto Min(tensor::type a, PDiff b) { return Min(Broadcast(Const(a), b->shape), b); }
inline auto Min(PDiff a, tensor::type b) { return Min(a, Broadcast(Const(b), a->shape)); }

struct MaxT : public Diff_vv {
    MaxT(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override { EACH(v) v[i] = max(va[i], vb[i]); }
    void Backward() override {
        EACH(ga) ga[i] += (va[i] > vb[i]) * g[i];
        EACH(gb) gb[i] += (vb[i] > va[i]) * g[i];
    }
};
Declare2(Max, MaxT);
inline auto Max(tensor::type a, PDiff b) { return Max(Broadcast(Const(a), b->shape), b); }
inline auto Max(PDiff a, tensor::type b) { return Max(a, Broadcast(Const(b), a->shape)); }

struct Concat : public Diff2 {
    Concat(PDiff a, PDiff b) : Diff2(a, b) {
        // TODO generalize for more dimensions
        Check(a->rank == 1);
        Check(b->rank == 1);
        Reshape({uint(a->size + b->size)});
    }

    void Forward() override {
        EACH(v) v[i] = va[i];
        EACH(v) v[i + va.size] = vb[i];
    }

    void Backward() override {
        EACH(ga) ga[i] += g[i];
        EACH(gb) gb[i] += g[i + va.size];
    }
};

// TODO Splice tensor operator

struct Conv1D : public Diff2 {
    Conv1D(PDiff a, PDiff b, int offset) : Diff2(a, b), offset(offset) {
        Check(a->rank == 1);
        Check(b->rank == 1);
        Reshape(a->shape);
    }

    void Forward() override {
        EACH(v) {
            tensor::type sum = 0;
            // TODO check
            size_t begin = (offset > i) ? offset - i : 0;
            size_t end = min<size_t>(vb.size, offset - i + va.size);
            for (size_t j = begin; j < end; j++) sum += va[i + j - offset] * vb[j];
            v[i] = sum;
        }
    }

    void Backward() override { Check(false, "not implemented"); }

    int offset;
};

struct Deconv1D : public Diff2 {
    Deconv1D(PDiff a, PDiff b) : Diff2(a, b) {}
};

struct MaxPool1D : public Diff1 {
    MaxPool1D(PDiff a) : Diff1(a) {
        Check(a->rank == 1);
        const uint m = (a->shape(0) + 1) / 2;
        Reshape({m});
    }

    void Forward() override {
        const size_t m = a->shape(0) / 2;
        for (size_t i = 0; i < m; i++) v[i] = max(va[i * 2], va[i * 2 + 1]);
        if (a->shape(0) % 2 == 1) v[m] = va[m * 2];
    }

    void Backward() override {
        if (!ga) return;
        Check(false);
        const size_t m = shape(0);
        for (size_t i = 0; i < m; i++) g[i] = max(va[i * 2], va[i * 2 + 1]);
    }
};

struct MaxPool2D : public Diff1 {
    MaxPool2D(PDiff a);
    void Forward() override;
    void Backward() override;
};

struct Reshape : public Diff1 {
    Reshape(PDiff a, tensor_shape shape) : Diff1(a) {
        v.reshape(shape);
        Check(a->size == v.size);
    }

    void Forward() override { EACH(v) v[i] = va[i]; }
    void Backward() override { EACH(ga) ga[i] = g[i]; }
};

// m can be zero, one or more dimensions
// p and q are exactly one dimension
// v{m,p} x m{q,p} -> {m,q}
struct VecMatMulT : public Diff2 {
    VecMatMulT(PDiff a, PDiff b) : Diff2(a, b) {
        tensor_shape as = a->shape, bs = b->shape;
        Check(as.size > 0);
        Check(bs.size == 2);
        Check(as.back() == bs.back());
        Reshape(as.pop_back().push_back(bs[0]));
    }

    void Forward() override {
        const size_t m = va.shape().pop_back().volume;
        const size_t q = vb.shape()[0];
        const size_t p = vb.shape()[1];

        for (auto im : range(m)) {
            for (auto iq : range(q)) {
                tensor::type s = 0;
                for (auto ip : range(p)) s += va[im * p + ip] * vb[iq * p + ip];
                v[im * q + iq] = s;
            }
        }
    }

    void Backward() override {
        const size_t m = va.shape().pop_back().volume;
        const size_t q = vb.shape()[0];
        const size_t p = vb.shape()[1];

        if (ga) {
            for (auto im : range(m))
                for (auto ip : range(p)) {
                    tensor::type s = ga[im * p + ip];
                    for (auto iq : range(q)) s += g[im * q + iq] * vb[iq * p + ip];
                    ga[im * p + ip] = s;
                }
        }
        if (gb) {
            for (auto iq : range(q))
                for (auto ip : range(p)) {
                    tensor::type s = gb[iq * p + ip];
                    for (auto im : range(m)) s += g[im * q + iq] * va[im * p + ip];
                    gb[iq * p + ip] = s;
                }
        }
    }
};

Declare2(VecMatMul, VecMatMulT);

struct SumT : public Diff1 {
    SumT(PDiff a) : Diff1(a) { Reshape({1}); }
    void Forward() override { v[0] = Sum(va); }
    void Backward() override { EACH(ga) ga[i] += g[0]; }
};

Declare1(Sum);

// TODO Unary horizontal Min and Max!

struct RunningAverageT : public Diff1 {
    RunningAverageT(PDiff a, float k) : Diff1(a), k(k) { v.reshape(a->shape); }
    void Forward() override { EACH(v) v[i] = v[i] * (1 - k) + va[i] * k; }
    tensor::type k;
};

inline PDiff RunningAverage(PDiff a, float k) { return make_shared<RunningAverageT>(a, k); }

// Returns mean of all input values from prev epoch OR 0 if first epoch.
struct EpochMeanT : public Diff1 {
    EpochMeanT(PDiff a) : Diff1(a) { v.reshape({1}); }
    void Forward() override {
        EACH(va) sum += va[i];
        count += va.size;
    }
    void EndEpoch() override {
        v[0] = sum / count;
        sum = 0;
        count = 0;
    }

    double sum = 0;
    uint count = 0;
};

inline PDiff EpochMean(PDiff a, float init) {
    auto p = make_shared<EpochMeanT>(a);
    p->v[0] = init;
    return p;
}

inline PDiff Mean(PDiff a) { return Sum(a) / a->size; }

inline PDiff MeanSquareError(PDiff a, PDiff b) { return Mean(Sqr(a - b)); }

template <typename T>
inline T RelativeError(T a, T b) {
    auto abs_error = abs(a - b);
    return abs_error / max(abs(a), abs(b));
}

template <typename T>
inline bool EqualFP(T a, T b) {
    auto abs_error = abs(a - b);
    if (a == 0 || b == 0) return abs_error <= 1e-6;
    return abs_error / max(abs(a), abs(b)) <= 1e-6;
}

template <typename T>
inline void CheckEqualFP(T a, T b) {
    if (!EqualFP(a, b))
        Fail(format("Expected equality of %.10f and %.10f (abs error %g, rel error %g)", a, b, abs(a - b),
                    RelativeError(a, b)));
}

struct ValueCmpT : public Diff_vv {
    ValueCmpT(PDiff a, PDiff b) : Diff_vv(a, b) {}
    void Forward() override {
        EACH(v) {
            v[i] = va[i];
            CheckEqualFP(va[i], vb[i]);
        }
    }
    void Backward() override {
        EACH(ga) ga[i] = g[i];
        EACH(gb) gb[i] = g[i];
    }
};

Declare2(ValueCmp, ValueCmpT);

struct GradCmpT : public DiffA {
    GradCmpT(PDiff a) : DiffA(a) {}
    void Forward() override { EACH(v) v[i] = va[i]; }
    void Backward() override {
        if (++counter != other->counter) return;
        EACH(ga) ga[i] = g[i];
        EACH(g) CheckEqualFP(g[i], other->g[i]);
    }
    void EndEpoch() override {
        if (counter != other->counter) Fail(format("%s != %s", counter, other->counter));
    }

    uint counter = 0;
    GradCmpT* other;  // can't use shared_ptr because cyclic references
};

inline pair<PDiff, PDiff> GradCmp(PDiff a) {
    auto p = make_shared<GradCmpT>(a);
    auto q = make_shared<GradCmpT>(a);
    p->other = q.get();
    q->other = p.get();
    return {p, q};
}

struct BinaryCrossEntropyT : public Diff2 {
    constexpr static tensor::type eps = 1e-6, one = 1;

    BinaryCrossEntropyT(PDiff ref, PDiff out) : Diff2(ref, out) {
        Check(a->shape == b->shape);
        Check(out->g);
        Check(!ref->g);
        Reshape({1});
    }
    void Forward() override {
        double s = 0;
        EACH(va) {
            s -= va[i] * std::log(max(eps, vb[i]));
            s -= (one - va[i]) * std::log(max(eps, one - vb[i]));
        }
        v[0] = s / a->size;
    }
    void Backward() override {
        Check(g[0] == 1);
        EACH(gb) {
            // TODO try removing the max
            auto p = va[i] / max(eps, vb[i]) * (eps < vb[i]);
            auto q = (one - va[i]) / max(eps, one - vb[i]) * (eps < one - vb[i]);
            gb[i] -= (p - q) / a->size; // not multiplying with g[0] == 1
        }
    }
};

template<bool Checker = false>
inline PDiff BinaryCrossEntropy(PDiff ref, PDiff out) {
    if (Checker) {
        auto [out0, out1] = GradCmp(out);
        auto one = Mean(-(ref * Log(Max(1e-6, out0)) + (1 - ref) * Log(Max(1e-6, 1 - out0))));
        auto two = make_shared<BinaryCrossEntropyT>(ref, out1);
        return ValueCmp(one, two);
    }
    return make_shared<BinaryCrossEntropyT>(ref, out);
}

struct AccuracyT : public Diff2 {
    AccuracyT(PDiff ref, PDiff out) : Diff2(ref, out) {
        Check(a->shape == b->shape);
        Check(!ref->g);
        v.reshape({1});
    }
    void Forward() override {
        double s = 0;
        EACH(va) s += std::abs(va[i] - vb[i]) < 0.5f;
        v[0] = s / a->size;
    }
};

Declare2(Accuracy, AccuracyT);

inline PDiff Softmax(PDiff a) {
    PDiff e = Exp(a);
    return e / Sum(e) << "softmax";
}

inline PDiff BatchNorm(PDiff a, tensor::type k) {
    auto mean = Mean(a) << "bn_mean";
    auto stdev = Sqrt(Mean(Sqr(a - mean))) + k << "bn_stdev";

    auto scale = 1 / stdev << "bn_scale";
    auto shift = mean / stdev << "bn_shift";
    auto norm = a * scale - shift << "bn_norm";

    auto gamma = Param(a->shape, 1) << "bn_gamma";
    auto beta = Param(a->shape, 0) << "bn_beta";
    return norm * gamma + beta << "bn";
}

inline PDiff FullyConnectedNB(PDiff a, uint size, shared_ptr<Init> w_init = make_shared<Init>()) {
    auto w = Param({size, a->shape->back()}, w_init) << "fc_w";
    return VecMatMul(a, w) << "fc";
}

inline PDiff FullyConnected(PDiff a, uint size, shared_ptr<Init> w_init = make_shared<Init>()) {
    auto w = Param({size, a->shape->back()}, w_init) << "fc_w";
    auto b = Param({size}) << "fc_b";
    return VecMatMul(a, w) + b << "fc";
}
