#pragma once
#include <core/diff.h>
#include <core/std.h>

struct ValueFunction : public Model {
    ValueFunction(uint16_t input_size) {
        input = Data({input_size}, "input");
        reference = Data({1}, "reference");
        NormalInit w_init(0.01f);
        auto fc1 = Relu(FullyConnected(input, 128, w_init));
        output = FullyConnected(fc1, 1, w_init);
        loss = MeanSquareError(output, reference);
    }

    float Loss(const tensor in, const tensor ref) {
        input->v = in; // TODO don't allow changing chape!
        reference->v = ref;
        Forward();
        return loss->v[0];
    }

    float Predict(const tensor in) {
        input->v = in;
        Forward(); // TODO no need to compute loss here
        return output->v[0];
    }
};
