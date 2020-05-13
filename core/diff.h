#pragma once
#include <core/std.h>
#include <core/tensor.h>

struct Diff {
    Diff(tensor_shape shape) : m_y(shape), m_dy(shape) {}

    const tensor y() const { return m_y; }
    tensor dy() { return m_dy; }

    auto shape() const { return m_y.shape(); }
    size_t size() const { return m_y.size(); }

    void Forward() {}
    void Backward() {}

   protected:
    vtensor m_y, m_dy;
};

struct Constant : public Diff {};

// learnable parameter, can be saved
struct Parameter : public Diff {};

// input and reference value
struct Variable : public Diff {
    Variable(tensor_shape shape) : Diff(shape) {}

    void set(const tensor y) {
        Check(m_y.shape() == y.shape());
        m_y = y;
    }
};

struct DiffFunc : public Diff {
    DiffFunc(Diff* input) : Diff(input->shape()), m_input(input) {}

    void Forward() { Func(m_input->y(), m_y); }

    void Backward() {
        Grad(m_input->y(), m_y, m_input->dy());
        for (size_t i = 0; i < m_dy.size(); i++) m_input->dy()[i] *= m_dy[i];
    }

    virtual void Func(const tensor x, tensor y) = 0;
    virtual void Grad(const tensor x, const tensor y, tensor dy) = 0;

   private:
    Diff* m_input;
};

struct Relu : public DiffFunc {
    Relu(Diff* input) : DiffFunc(input) {}

    void Func(const tensor x, tensor y) override {
        for (size_t i = 0; i < y.size(); i++) y[i] = (x[i] > 0) ? x[i] : 0;
    }

    void Grad(const tensor x, const tensor y, tensor dy) override {
        for (size_t i = 0; i < y.size(); i++) dy[i] = (x[i] > 0) ? 1 : 0;
    }
};

struct Sigmoid : public DiffFunc {
    Sigmoid(Diff* input) : DiffFunc(input) {}

    void Func(const tensor x, tensor y) override {
        for (size_t i = 0; i < y.size(); i++) y[i] = 1 / (1 + exp(-x[i]));
    }

    void Grad(const tensor x, const tensor y, tensor dy) override {
        for (size_t i = 0; i < y.size(); i++) dy[i] = y[i] * (1 - y[i]);
    }
};

inline tensor::type Dot(const tensor a, const tensor b) {
    Check(a.shape() == b.shape());
    tensor::type sum = 0;
    for (size_t i = 0; i < a.size(); i++) sum += a[i] * b[i];
    return sum;
}

struct FullyConnected : public Diff {
    FullyConnected(Diff* input, uint16_t size, float variance, std::mt19937& random)
        : Diff(Shape(size)), m_input(input), m_w(Concat(size, input->shape())), m_b(Shape(size), 0) {
        std::normal_distribution<float> dis(0.0f, variance);
        for (size_t i = 0; i < m_w.size(); i++) m_w[i] = dis(random);
    }

    void Forward() {
        const tensor x = m_input->y();
        const size_t m = x.size();
        const size_t n = m_b.size();

        for (size_t i = 0; i < n; i++) {
            float s = m_b[i];
            // TODO vectorize
            for (size_t j = 0; j < m; j++) s += m_w[i * m + j] * x[j];
            m_y[i] = s;
        }
    }

    void Backward() {
        const tensor x = m_input->y();
        tensor dx = m_input->dy();
        const size_t m = x.size();
        const size_t n = m_b.size();

        for (size_t j = 0; j < m; j++) {
            float s = 0;
            for (size_t i = 0; i < n; i++) s += m_dy[i] * m_w[i * m + j];
            dx[j] = s;
        }
    }

   private:
    Diff* m_input;
    vtensor m_w, m_b;
};

struct MeanSquareError : public Diff {
    MeanSquareError(Diff* input, const Diff* reference) : Diff({1}), m_input(input), m_reference(reference) {
        Check(input->shape() == reference->shape());
    }

    void Forward() {
        const tensor x = m_input->y();
        const tensor r = m_reference->y();
        float error = 0;
        for (size_t i = 0; i < x.size(); i++) error += sqr(x[i] - r[i]);
        m_y[0] = error;
    }

    void Backward() {
        const tensor x = m_input->y();
        const tensor r = m_reference->y();
        tensor dx = m_input->dy();

        for (size_t i = 0; i < x.size(); i++) dx[i] = 2 * (x[i] - r[i]);
    }

   private:
    const Diff* m_reference;
    Diff* m_input;
};

struct BatchNorm : public Diff {
    BatchNorm(Diff* input, float delta) : Diff(input->shape()), m_delta(delta) {}

    void Forward() {
        const auto& x = m_input->y();

        const float mean = Sum<float>(x) / x.size();

        float s = 0;
        for (size_t i = 0; i < x.size(); i++) s += sqr(x[i] - mean);
        const float stdev = sqrt(m_delta + s / x.size());

        for (size_t i = 0; i < x.size(); i++) m_y[i] = (x[i] - mean) / stdev;
    }

    void Backward() {
        const tensor x = m_input->y();
        tensor dx = m_input->dy();

        // TODO(Marko)
        Check(false, "unfinished");
        // dj/dx  = dj/dy * d (1/stdev )/dx
    }

   private:
    const float m_delta;
    Diff* m_input;
};

template <typename Model>
void TrainWithSGD(Model& model, const tensor in, tensor out, std::mt19937& random) {
    Check(in.shape().size() > 0);
    Check(out.shape().size() > 0);
    Check(in.shape().back() == out.shape().back());

    Check(PopBack(in.shape()) == model.input.shape());
    Check(PopBack(out.shape()) == model.output.shape());

    vector<uint32_t> samples;
    samples.resize(in.shape().back());
    for (uint32_t i = 0; i < samples.size(); i++) samples[i] = i;
    std::shuffle(samples.begin(), samples.end(), random);

    for (uint32_t index : samples) {
        // model.Train(in.sub(index), out.sub(index));
        // x
        // y = f(x, w)
        // loss = norm(y - ref)
        // optimize loss == 0  ==>  norm(f(x, w) - ref) == 0  ==>  w -= alpha *  d/dw f(x, w)
        // w -= alpha * d(x, w)
    }
}
