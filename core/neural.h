#pragma once
#include <core/std.h>
#include <core/tensor.h>

class Layer {
   public:
    Layer(cspan<uint32_t> shape) : m_y(shape), m_dy(shape) {}

    const Tensor<float>& y() const { return m_y; }
    Tensor<float>& dy() { return m_dy; }

    auto shape() const { return m_y.shape(); }
    size_t size() const { return m_y.size(); }

    void Forward() {}
    void Backward() {}

   protected:
    Tensor<float> m_y, m_dy;
};

class InputLayer : public Layer {
   public:
    InputLayer(cspan<uint32_t> shape) : Layer(shape) {}
    void set(const TensorSpan<const float>& y) {
        Check(m_y.shape() == y.shape());
        m_y = y;
    }
};

class ReluLayer : public Layer {
   public:
    ReluLayer(Layer* input) : Layer(input->shape()), m_input(input) {}

    void Forward() {
        const auto& x = m_input->y();
        for (size_t i = 0; i < m_y.size(); i++) m_y[i] = (x[i] > 0) ? x[i] : 0;
    }

    void Backward() {
        const auto& x = m_input->y();
        auto& dx = m_input->dy();
        for (size_t i = 0; i < x.size(); i++) dx[i] = (x[i] > 0) ? m_dy[i] : 0.0f;
    }

   private:
    Layer* m_input;
};

inline double Sigmoid(double x) { return 1 / (1 + exp(-x)); }

class SigmoidLayer : public Layer {
   public:
    SigmoidLayer(Layer* input) : Layer(input->shape()), m_input(input) {}

    void Forward() {
        const auto& x = m_input->y();
        for (size_t i = 0; i < x.size(); i++) m_y[i] = Sigmoid(x[i]);
    }

    void Backward() {
        const auto& x = m_input->y();
        auto& dx = m_input->dy();
        for (size_t i = 0; i < x.size(); i++) dx[i] = m_dy[i] * m_y[i] * (1 - m_y[i]);
    }

   private:
    Layer* m_input;
};

inline float Dot(const TensorSpan<const float>& a, const TensorSpan<const float>& b) {
    Check(a.shape() == b.shape());
    float sum = 0;
    for (size_t i = 0; i < a.size(); i++) sum += a[i] * b[i];
    return sum;
}

vector<uint32_t> Concat(uint32_t a, cspan<uint32_t> b) {
    vector<uint32_t> out(1 + b.size());
    out[0] = a;
    Copy(b, span<uint32_t>(out.data() + 1, out.size() - 1));
    return out;
}

class FullyConnectedLayer : public Layer {
   public:
    FullyConnectedLayer(Layer* input, uint32_t size, float variance, std::mt19937& random)
        : Layer({size}), m_input(input), m_w(Concat(size, input->shape())), m_b({size}, 0.0f) {
        std::normal_distribution<float> dis(0.0f, variance);
        for (size_t i = 0; i < m_w.size(); i++) m_w[i] = dis(random);
    }

    void Forward() {
        const auto& x = m_input->y();
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
        const auto& x = m_input->y();
        auto& dx = m_input->dy();
        const size_t m = x.size();
        const size_t n = m_b.size();

        for (size_t j = 0; j < m; j++) {
            float s = 0;
            for (size_t i = 0; i < n; i++) s += m_dy[i] * m_w[i * m + j];
            dx[j] = s;
        }
    }

   private:
    Layer* m_input;
    Tensor<float> m_w, m_b;
};

class MeanSquareErrorLayer : public Layer {
   public:
    MeanSquareErrorLayer(Layer* input, const Layer* reference) : Layer({1}), m_input(input), m_reference(reference) {
        Check(input->shape() == reference->shape());
    }

    void Forward() {
        const auto& x = m_input->y();
        const auto& r = m_reference->y();
        float error = 0;
        for (size_t i = 0; i < x.size(); i++) {
            float e = x[i] - r[i];
            error += e * e;
        }
        m_y[0] = error;
    }

    void Backward() {
        const auto& x = m_input->y();
        const auto& r = m_reference->y();
        auto& dx = m_input->dy();

        for (size_t i = 0; i < x.size(); i++) dx[i] = 2 * (x[i] - r[i]);
    }

   private:
    const Layer* m_reference;
    Layer* m_input;
};

template <typename Model>
void TrainWithSGD(Model& model, const Tensor<float>& in, const Tensor<float>& out, std::mt19937& random) {
    Check(in.shape().size() > 0);
    Check(out.shape().size() > 0);
    Check(in.shape().back() == out.shape().back());

    Check(in.shape().pop_back() == model.input.shape());
    Check(out.shape().pop_back() == model.output.shape());

    vector<uint32_t> samples;
    samples.resize(in.shape().back());
    for (uint32_t i = 0; i < samples.size(); i++) samples[i] = i;
    std::shuffle(samples.begin(), samples.end(), random);

    for (uint32_t index : samples) {
        model.Train(in.sub(index), out.sub(index));
        // x
        // y = f(x, w)
        // loss = norm(y - ref)
        // optimize loss == 0  ==>  norm(f(x, w) - ref) == 0  ==>  w -= alpha *  d/dw f(x, w)
        // w -= alpha * d(x, w)
    }
}
