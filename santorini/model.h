#pragma once
#include <core/std.h>
#include <core/neural.h>

struct FeedForwardNetwork {
    InputLayer input;
    InputLayer reference;

    // hidden
    FullyConnectedLayer fc1;
    ReluLayer relu;
    FullyConnectedLayer fc2;

    SigmoidLayer output;
    MeanSquareErrorLayer loss;

    FeedForwardNetwork(uint32_t input_size, uint32_t fc1_size, std::mt19937& random)
        : input({input_size}),
          reference({1}),
          fc1(&input, fc1_size, 0.01f, random),
          relu(&fc1),
          fc2(&relu, 1, 0.01f, random),
          output(&fc2),
          loss(&output, &reference) {}

    void Load(istream& is) {
        // is >> input >> fc1 >> relu >> fc2 >> sigmoid >> loss;
    }

    void Save(ostream& os) const {
        // os << input << fc1 << relu << fc2 << sigmoid << loss;
    }

    void Forward() {
        fc1.Forward();
        relu.Forward();
        fc2.Forward();
        output.Forward();
    }

    void Backward() {
        output.Backward();
        fc2.Backward();
        relu.Backward();
        fc1.Backward();
    }

    void Train(const TensorSpan<const float>& in, const TensorSpan<const float>& ref) {
        input.set(in);
        reference.set(ref);

        Forward();
        loss.Forward();
        // TODO how to propagate loss back in network?
        loss.Backward();
        Backward();
    }

    void TrainBatch(const TensorSpan<const float>& in, const TensorSpan<const float>& ref) {
        // TODO(Marko)
    }

    float Loss(const TensorSpan<const float>& in, const TensorSpan<const float>& ref) {
        input.set(in);
        reference.set(ref);
        Forward();
        loss.Forward();
        return loss.y()[0];
    }

    float Predict(const Tensor<float>& in) {
        input.set(in);
        Forward();
        return output.y()[0];
    }
};
