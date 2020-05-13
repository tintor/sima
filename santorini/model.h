#pragma once
#include <core/neural.h>
#include <core/std.h>

struct FeedForwardNetwork {
    Input input;
    Input reference;

    // hidden
    FullyConnected fc1;
    Relu relu;
    FullyConnected fc2;

    Sigmoid output;
    MeanSquareError loss;

    FeedForwardNetwork(uint16_t input_size, uint16_t fc1_size, std::mt19937& random)
        : input({input_size, 0, 0, 0}),
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

    void Train(const tensor in, const tensor ref) {
        input.set(in);
        reference.set(ref);

        Forward();
        loss.Forward();
        // TODO how to propagate loss back in network?
        loss.Backward();
        Backward();
    }

    void TrainBatch(const tensor in, const tensor ref) {
        // TODO(Marko)
    }

    float Loss(const tensor in, const tensor ref) {
        input.set(in);
        reference.set(ref);
        Forward();
        loss.Forward();
        return loss.y()[0];
    }

    float Predict(const tensor in) {
        input.set(in);
        Forward();
        return output.y()[0];
    }
};
