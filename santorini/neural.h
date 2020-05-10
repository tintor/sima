#pragma once
#include <core/std.h>
#include <santorini/tensor.h>

class Layer {
   public:
    Layer(cspan<uint32_t> shape) : m_output(shape) {}

    const Tensor<float>& Output() const { return m_output; }
    size_t Size() const { return m_output.size(); }
    float operator[](size_t index) const { return m_output[index]; }
    float operator()(cspan<uint32_t> index) const { return m_output(index); }

    void Forward() {}
    void Backward() {}

   protected:
    Tensor<float> m_output;
};

class InputLayer : public Layer {
   public:
    InputLayer(cspan<uint32_t> shape) : Layer(shape) {}
    float& operator()(cspan<uint32_t> index) { return m_output(index); }
    auto& Tensor() { return m_output; }
};

inline double Relu(double a) { return a > 0 ? a : 0; }

class ReluLayer : public Layer {
   public:
    ReluLayer(const Layer* input) : Layer(input->Output().shape()), m_input(input) {}

    void Forward() {
        for (size_t i = 0; i < m_output.size(); i++) m_output[i] = Relu(m_input->Output()[i]);
    }

    void Backward() {}

   private:
    const Layer* m_input;
};

inline double Sigmoid(double a) { return 1 / (1 + exp(a)); }

class SigmoidLayer : public Layer {
   public:
    SigmoidLayer(const Layer* input) : Layer(input->Output().shape()), m_input(input) {}

    void Forward() {
        const auto& input = m_input->Output();
        for (size_t i = 0; i < input.size(); i++) m_output[i] = Sigmoid(input[i]);
    }

    void Backward() {}

   private:
    const Layer* m_input;
};

class FullyConnectedLayer : public Layer {
   public:
    FullyConnectedLayer(const Layer* input, uint32_t size, float variance, std::mt19937& random)
        : Layer({size}), m_input(input), m_weights({uint32_t(size * input->Output().size())}), m_biases({size}) {
        std::normal_distribution<float> dis(0.0f, variance);
        for (size_t i = 0; i < m_weights.size(); i++) m_weights[i] = dis(random);
        for (size_t i = 0; i < m_biases.size(); i++) m_biases[i] = 0.0f;
    }

    void Forward() {
        const auto& input = m_input->Output();
        const size_t m = input.size();
        const size_t n = m_biases.size();

        for (size_t i = 0; i < n; i++) {
            float v = m_biases[i];
            // TODO vectorize
            for (size_t j = 0; j < m; j++) v += m_weights[i * m + j] * input[j];
            m_output[i] = v;
        }
    }

    void Backward() {}

   private:
    const Layer* m_input;
    Tensor<float> m_weights, m_biases;
};

struct FeedForwardNetwork {
    InputLayer input;
    FullyConnectedLayer fc1;
    ReluLayer relu;
    FullyConnectedLayer fc2;
    SigmoidLayer sigmoid;

    FeedForwardNetwork(uint32_t input_size, uint32_t fc1_size, std::mt19937& random)
        : input({input_size}),
          fc1(&input, fc1_size, 0.01f, random),
          relu(&fc1),
          fc2(&relu, 1, 0.01f, random),
          sigmoid(&fc2) {}

    void Load(istream& is) {
        // is >> input >> fc1 >> relu >> fc2 >> sigmoid;
    }

    void Save(ostream& os) const {
        // os << input << fc1 << relu << fc2 << sigmoid;
    }

    void Forward() {
        input.Forward();
        fc1.Forward();
        relu.Forward();
        fc2.Forward();
        sigmoid.Forward();
    }

    void Backward() {
        sigmoid.Backward();
        fc2.Backward();
        relu.Backward();
        fc1.Backward();
        input.Backward();
    }

    void TrainBatch1(const Tensor<float>& in, const Tensor<float>& out) {
        Check(out.shape().size() == 1);
        Check(out.shape()[0] == 1);
        input.Tensor() = in;
        Forward();
        float loss = out[0] - sigmoid[0];
        // TODO how to propagate loss back in network?
        Backward();
    }

    void TrainBatchN(const Tensor<float>& in, const Tensor<float>& out) {}

    float Predict(const Tensor<float>& in) {
        input.Tensor() = in;
        Forward();
        return sigmoid[0];
    }
};
