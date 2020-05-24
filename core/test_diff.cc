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

void Minimize(PDiff loss, float alpha, size_t iterations) {
    auto nodes = TopoSort({loss});
    Iterate(nodes, alpha, iterations);
}

TEST_CASE("diff: minimize circle", "[diff]") {
    auto x = Param({1}) << "x";
    auto y = Param({1}) << "y";
    x->v[0] = 3.14f;
    y->v[0] = 2.16f;
    Minimize(Sqr(x) + Sqr(y), 0.1, 70);

    REQUIRE(abs(x->v[0]) < 1e-6);
    REQUIRE(abs(y->v[0]) < 1e-6);
}

TEST_CASE("diff: minimize binary cross entropy ", "[diff]") {
    auto x = Param({1}) << "x";
    auto ref = Const(0);
    auto b = Logistic(x) << "b";
    //auto loss = BinaryCrossEntropy(ref, b);

    auto grad_cmp = GradCmp(b);
    auto loss0 = XBinaryCrossEntropy(ref, grad_cmp.first) << "loss0";
    auto loss1 = BinaryCrossEntropy(ref, grad_cmp.second) << "loss1";
    auto loss = ValueCmp(loss0, loss1);

    for (float e : range(-1.f, 1.f, 0.1f)) {
        x->v[0] = e;
        ref->v[0] = 0;
        auto nodes = TopoSort({loss});
        Forward(nodes);
        println("x:%s r:%s b:%s", x->v[0], ref->v[0], b->v[0]);
        Minimize(loss, 0.05, 500);
        REQUIRE(b->v[0] <= 0.05);
    }

    for (float e : range(-1.f, 1.f, 0.1f)) {
        x->v[0] = e;
        ref->v[0] = 1;
        auto nodes = TopoSort({loss});
        Forward(nodes);
        println("x:%s r:%s b:%s", x->v[0], ref->v[0], b->v[0]);
        Minimize(loss, 0.05, 500);
        REQUIRE(b->v[0] >= 0.95);
    }
}

TEST_CASE("diff: minimize booth", "[diff]") {
    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto x = Param({1}, init) << "x";
    auto y = Param({1}, init) << "y";
    Minimize(Sqr(x + 2 * y - 7) + Sqr(2 * x + y - 5), 0.01, 1000);

    REQUIRE(abs(x->v[0] - 1) <= 2e-5);
    REQUIRE(abs(y->v[0] - 3) <= 2e-5);
}

TEST_CASE("diff: minimize rastrigin", "[diff]") {
    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto x = Param({1}, init) << "x";
    auto y = Param({1}, init) << "y";
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

TEST_CASE("diff: learn perceptron, plane in 2d", "[diff_x]") {
    const int Batch = 20;
    auto x = Data({Batch, 1}) << "x";
    auto y = Data({Batch, 1}) << "y";
    auto ref = Data({Batch, 1}) << "ref";

    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto a = Param({1}, init) << "a";
    auto b = Param({1}, init) << "b";
    auto c = Param({1}, init) << "c";
    auto w = Param({1}, init) << "w";
    auto e = Param({1}, init) << "e";

    auto out = Logistic(x * a + y * b + c, 15) << "out";

    // Diffing!
    //auto loss = XBinaryCrossEntropy(ref, out) << "loss0";
    auto grad_cmp = GradCmp(out);
    auto loss0 = XBinaryCrossEntropy(ref, grad_cmp.first) << "loss0";
    auto loss1 = BinaryCrossEntropy(ref, grad_cmp.second) << "loss1";
    auto loss = ValueCmp(loss0, loss1);

    auto accuracy = Accuracy(ref, out);

    Model model(loss, accuracy);
    model.Print();

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
    for (auto i : range(100)) metrics = model.Epoch(dataset, 0.1, random, true, i);
    Print(model.nodes);

    //REQUIRE(metrics.at("accuracy") >= 0.9960);
}

PDiff Neuron(PDiff x, PDiff y, shared_ptr<Init> init) {
    auto a = Param({1}, init) << "a";
    auto b = Param({1}, init) << "b";
    auto c = Param({1}, init) << "c";
    return Logistic(x * a + y * b + c);
}

PDiff Neuron(PDiff x, PDiff y, PDiff z, shared_ptr<Init> init) {
    auto a = Param({1}, init) << "a";
    auto b = Param({1}, init) << "b";
    auto c = Param({1}, init) << "c";
    auto d = Param({1}, init) << "d";
    return Logistic(x * a + y * b + z * c + d);
}

TEST_CASE("diff: learn two layer network, circle in 2d", "[diff]") {
    const int Batch = 20;
    auto x = Data({Batch, 1}) << "x";
    auto y = Data({Batch, 1}) << "y";
    auto ref = Data({Batch, 1}) << "ref";

    std::mt19937_64 random(1);
    auto init = make_shared<NormalInit>(1, random);
    auto h1 = Neuron(x, y, init) << "h1";
    auto h2 = Neuron(x, y, init) << "h2";
    auto h3 = Neuron(x, y, init) << "h3";
    auto out = Neuron(h1, h2, h3, init) << "out";

    auto loss = MeanSquareError(out, ref) << "loss";
    auto accuracy = Mean(Abs(out - ref) < 0.5) << "accuracy";
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

TEST_CASE("diff: learn FC perceptron, plane in 3d", "[diff]") {
    const int Batch = 40;
    auto in = Data({Batch, 3}) << "in";
    auto ref = Data({Batch, 1}) << "ref";

    std::mt19937_64 random(3);
    auto init = make_shared<NormalInit>(1, random);
    auto fc = FullyConnected(in, 1, init);
    auto out = Tanh(fc) * 0.5 + 0.5 << "out";

    auto loss = MeanSquareError(out, ref) << "loss";
    auto accuracy = Mean(Abs(out - ref) < 0.5) << "accuracy";
    Model model(loss, accuracy);

    // dataset
    const float A = 0.2, B = 0.4, C = -0.8, D = 0.1;
    UniformInit gen(-1, 1, random);
    const int Classes = 2;
    const int SamplesPerClass = 20000;
    const int Samples = Classes * SamplesPerClass;
    vtensor data_in({Samples, 3}, 0);
    vtensor data_ref({Samples, 1}, 0);
    int count[2] = {0, 0};
    int index = 0;
    while (index < Samples) {
        float dx = gen.get();
        float dy = gen.get();
        float dz = gen.get();
        float dr = (dx * A + dy * B + dz * C + D >= 0) ? 1 : 0;
        int c = round(dr);
        if (count[c] >= SamplesPerClass) continue;
        count[c] += 1;
        data_in(index, 0) = dx;
        data_in(index, 1) = dy;
        data_in(index, 2) = dz;
        data_ref[index] = dr;
        index += 1;
    }
    vector<pair<PDiff, tensor>> dataset = {{in, data_in}, {ref, data_ref}};

    // train!
    println("train ...");
    Print(model.nodes);
    Metrics metrics;
    for (auto i : range(1000)) {
        metrics = model.Epoch(dataset, 0.01, random, i % 100 == 99);
    }
    Print(model.nodes);

    REQUIRE(metrics.at("accuracy") >= 0.9975);
}

TEST_CASE("diff: fully connected, two layer network, circle in 2d", "[.][diff_p]") {
    const int Batch = 20;
    auto in = Data({Batch, 2}) << "in";
    auto ref = Data({Batch, 1}) << "ref";

    std::mt19937_64 random(2);
    auto init = make_shared<NormalInit>(1, random);

    //auto proc = BatchNorm(in, 1e-10) << "proc";

    auto hidden = Logistic(FullyConnected(in, 3, init)) << "hidden";
    // hidden = BatchNorm(hidden, 1e-10);
    auto out = Logistic(FullyConnected(hidden, 1, init)) << "out";

    auto loss = BinaryCrossEntropy(ref, out) << "loss";
    auto accuracy = Mean(Abs(out - ref) < 0.5) << "accuracy";
    Model model(loss, accuracy);

    // dataset
    UniformInit gen(-1, 1, random);
    const int Classes = 2;
    const int SamplesPerClass = 10000;
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

    metrics = model.Epoch(dataset, 0.01, random, true);
    Print(model.nodes);
    return;

    for (auto i : range(10000)) {
        metrics = model.Epoch(dataset, 0.01, random, i % 25 == 24);
        if (i % 25 == 24) Print(model.nodes);
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
