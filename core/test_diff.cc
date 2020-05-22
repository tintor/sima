#include <core/diff.h>

#include <catch.hpp>

void Iterate(span<PDiff> nodes, float alpha, size_t iterations) {
    for (size_t i = 0; i < iterations; i++) {
        Forward(nodes);
        ResetGradients(nodes);
        nodes.back()->g[0] = 1;
        Backward(nodes);
        GradientDescent(nodes, alpha);
    }
}

TEST_CASE("diff: minimize circle", "[diff]") {
    auto x = Param({1}, "x");
    auto y = Param({1}, "y");
    auto loss = Sqr(x) + Sqr(y);
    loss->name = "loss";

    auto nodes = TopoSort({loss});
    x->v[0] = 3.14f;
    y->v[0] = 2.16f;
    Forward(nodes);
    REQUIRE(loss->v[0] == sqr(3.14f) + sqr(2.16f));

    Iterate(nodes, 0.1, 70);
    Print(nodes);

    REQUIRE(abs(x->v[0]) < 1e-6);
    REQUIRE(abs(y->v[0]) < 1e-6);
}

void Minimize(PDiff loss, float alpha, size_t iterations) {
    auto nodes = TopoSort({loss});
    Iterate(nodes, alpha, iterations);
    Print(nodes);
}

TEST_CASE("diff: minimize booth", "[diff]") {
    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto x = Param({1}, "x", init);
    auto y = Param({1}, "y", init);
    Minimize(Sqr(x + 2 * y - 7) + Sqr(2 * x + y - 5), 0.01, 1000);

    REQUIRE(abs(x->v[0] - 1) <= 2e-5);
    REQUIRE(abs(y->v[0] - 3) <= 2e-5);
}

TEST_CASE("diff: minimize rastrigin", "[diff]") {
    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto x = Param({1}, "x", init);
    auto y = Param({1}, "y", init);
    Minimize(20 + (x * x - 10 * Cos(2 * PI * x)) + (y * y - 10 * Cos(2 * PI * y)), 0.01, 2000);

    // REQUIRE(abs(x->v[0]) < 1e-5);
    // REQUIRE(abs(y->v[0]) < 1e-5);
}

/*TEST_CASE("diff: minimize rosenbrock", "[diff]") {
    auto init = make_shared<NormalInit>(1, 0);
    auto x = Param({1}, "x", init);
    auto y = Param({1}, "y", init);
    Minimize(100*Sqr(y - x*x) + Sqr(1 - x), 0.001, 2000);

    REQUIRE(abs(x->v[0] - 1) < 1e-5);
    REQUIRE(abs(y->v[0] - 1) < 1e-5);
}*/

#if 0
    Optimize("easom", 100, return -cos(x)*cos(y)*exp(-sqr(x - PI) - sqr(y - PI)));
    Optimize("himmelblau", 5, return sqr(x*x + y - 11) + sqr(x + y*y - 7));
    Optimize("styblinski–tang", 5, return x*x*x*x - 16*x*x + 5*x + y*y*y*y - 16*y*y + 5*y);
    Optimize("schaffer n.2", 100, return 0.5 + (sqr(sin(x*x - y*y)) - 0.5) / sqr(1 + 0.001*(x*x + y*y)));
    Optimize("hölder table", 10, return -abs(sin(x)*cos(y)*exp(abs(1 - sqrt(x*x + y*y)/PI))));
#endif

TEST_CASE("diff: learn perceptron, plane in 2d", "[diff]") {
    const int Batch = 20;
    auto x = Data({Batch, 1}, "x");
    auto y = Data({Batch, 1}, "y");
    auto ref = Data({Batch, 1}, "ref");

    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto a = Param({1}, "a", init);
    auto b = Param({1}, "b", init);
    auto c = Param({1}, "c", init);
    auto w = Param({1}, "w", init);
    auto e = Param({1}, "e", init);

    auto out = Tanh(x * a + y * b + c) * 0.5 + 0.5;
    out->name = "out";
    auto loss = MeanSquareError(out, ref);
    loss->name = "loss";
    auto accuracy = Mean(Abs(out - ref) < 0.5);
    accuracy->name = "accuracy";
    Model model(loss, accuracy);

    // dataset
    const float A = 0.4, B = 0.6, C = -0.4;
    UniformInit gen(-1, 1, random);
    const int Classes = 2;
    const int SamplesPerClass = 20000;
    const int Samples = Classes * SamplesPerClass;
    vtensor data_x({Samples, 1}, 0);
    vtensor data_y({Samples, 1}, 0);
    vtensor data_r({Samples, 1}, 0);
    int count[2] = {0, 0};
    int index = 0;
    while (index < Samples) {
        float dx = gen.get();
        float dy = gen.get();
        float dr = (dx * A + dy * B + C >= 0) ? 1 : 0;
        int c = round(dr);
        if (count[c] >= SamplesPerClass) continue;
        count[c] += 1;
        data_x[index] = dx;
        data_y[index] = dy;
        data_r[index] = dr;
        index += 1;
    }
    vector<pair<PDiff, tensor>> dataset = {{x, data_x}, {y, data_y}, {ref, data_r}};

    // train!
    println("train ...");
    Print(model.nodes);
    Metrics metrics;
    for (auto i : range(800)) metrics = model.Epoch(dataset, 0.1, random, false);
    Print(model.nodes);

    REQUIRE(metrics.at("accuracy") >= 0.9965);
}

PDiff Neuron(PDiff x, PDiff y, string_view name, shared_ptr<Init> init) {
    auto a = Param({1}, "a", init);
    auto b = Param({1}, "b", init);
    auto c = Param({1}, "c", init);
    auto h = Logistic(x * a + y * b + c);
    h->name = name;
    return h;
}

PDiff Neuron(PDiff x, PDiff y, PDiff z, string_view name, shared_ptr<Init> init) {
    auto a = Param({1}, "a", init);
    auto b = Param({1}, "b", init);
    auto c = Param({1}, "c", init);
    auto d = Param({1}, "d", init);
    auto h = Logistic(x * a + y * b + z * c + d);
    h->name = name;
    return h;
}

TEST_CASE("diff: learn two layer network, circle in 2d", "[diff_circle]") {
    const int Batch = 20;
    auto x = Data({Batch, 1}, "x");
    auto y = Data({Batch, 1}, "y");
    auto ref = Data({Batch, 1}, "ref");

    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto h1 = Neuron(x, y, "h1", init);
    auto h2 = Neuron(x, y, "h2", init);
    auto h3 = Neuron(x, y, "h3", init);
    auto out = Neuron(h1, h2, h3, "out", init);

    auto loss = MeanSquareError(out, ref);
    loss->name = "loss";
    auto accuracy = Mean(Abs(out - ref) < 0.5);
    accuracy->name = "accuracy";
    Model model(loss, accuracy);

    // dataset
    UniformInit gen(-1, 1, random);
    const int Classes = 2;
    const int SamplesPerClass = 20000;
    const int Samples = Classes * SamplesPerClass;
    vtensor data_x({Samples, 1}, 0);
    vtensor data_y({Samples, 1}, 0);
    vtensor data_r({Samples, 1}, 0);
    int count[2] = {0, 0};
    int index = 0;
    while (index < Samples) {
        float dx = gen.get();
        float dy = gen.get();
        float dr = (dx * dx + dy * dy >= 0.5) ? 1 : 0;
        int c = round(dr);
        if (count[c] >= SamplesPerClass) continue;
        count[c] += 1;
        data_x[index] = dx;
        data_y[index] = dy;
        data_r[index] = dr;
        index += 1;
    }
    vector<pair<PDiff, tensor>> dataset = {{x, data_x}, {y, data_y}, {ref, data_r}};

    // train!
    println("train ...");
    Print(model.nodes);
    Metrics metrics;
    for (auto i : range(1000)) metrics = model.Epoch(dataset, 0.1, random, i % 25 == 24);
    Print(model.nodes);

    REQUIRE(metrics.at("accuracy") >= 0.9919);
}

TEST_CASE("diff: fully connected, two layer network, circle in 2d", "[diff_fc]") {
    const int Batch = 100;
    auto in = Data({Batch, 2}, "in");
    auto ref = Data({Batch, 1}, "ref");

    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto hidden = Logistic(FullyConnected(in, 2, init));
    auto out = Logistic(FullyConnected(hidden, 1, init));

    auto loss = MeanSquareError(out, ref);
    loss->name = "loss";
    auto accuracy = Mean(Abs(out - ref) < 0.5);
    accuracy->name = "accuracy";
    Model model(loss, accuracy);

    // dataset
    UniformInit gen(-1, 1, random);
    const int Classes = 2;
    const int SamplesPerClass = 20000;
    const int Samples = Classes * SamplesPerClass;
    vtensor data_in({Samples, 2}, 0);
    vtensor data_ref({Samples, 1}, 0);
    int count[2] = {0, 0};
    int index = 0;
    while (index < Samples) {
        float dx = gen.get();
        float dy = gen.get();
        float dr = (dx * dx + dy * dy >= 0.5) ? 1 : 0;
        int c = round(dr);
        if (count[c] >= SamplesPerClass) continue;
        count[c] += 1;
        data_in(index, 0) = dx;
        data_in(index, 1) = dy;
        data_ref[index] = dr;
        index += 1;
    }
    vector<pair<PDiff, tensor>> dataset = {{in, data_in}, {ref, data_ref}};

    // train!
    println("train ...");
    Print(model.nodes);
    Metrics metrics;
    for (auto i : range(1000)) {
        metrics = model.Epoch(dataset, 0.03, random, i % 25 == 24);
        if (!Bounded(model.nodes, 1e6)) break;
    }
    Print(model.nodes);

    REQUIRE(metrics.at("accuracy") >= 0.9919);
}

// Classification:
// learn hyper plane in Nd
// learn circle in 2d
// learn hyper sphere in Nd
// spiral in 2d

// Regression:
// - multiply two inputs: x * y
// - sqr
// - sqrt
// - sine in 1d
// - sine in 2d : sin(sqrt(x^2 + y^2)

// Neuro-plasticity:
// 1) train same model on datasets A and B separately, then take model trained of A and train it on B, should be able to
// get to similar performance level as when only trained on B 2) train fully on A, then fully on B, then some refresher
// on A, and make sure that both skills are retained

// Be able to detect if input is 1) in-domain (ie. similar to training data) or 2) out-of-domain (very different from
// training data)

// Idea: add 2nd order derivatives!
// Idea: solve motion planning problems! One layer/diff for one time step. All costs add up to loss. Parameters are
// control inputs! Idea: could take advantage of mini batches to solve multiple problems at once. Idea: would need to be
// able to compute Lp-distance on GPU

// Idea: let model decide if it wants to propagate gradients to any specific part of the network
