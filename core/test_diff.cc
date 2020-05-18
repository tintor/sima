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
    Setup(nodes);

    x->v[0] = 3.14f;
    y->v[0] = 2.16f;
    Forward(nodes);
    REQUIRE(loss->v[0] == sqr(3.14f) + sqr(2.16f));

    Iterate(nodes, 0.1, 70);
    Print(nodes, true, true);

    REQUIRE(abs(x->v[0]) < 1e-6);
    REQUIRE(abs(y->v[0]) < 1e-6);
}

void Minimize(PDiff loss, float alpha, size_t iterations) {
    auto nodes = TopoSort({loss});
    Setup(nodes);
    Iterate(nodes, alpha, iterations);
    Print(nodes, true, true);
}

TEST_CASE("diff: minimize booth", "[diff]") {
    auto init = make_shared<NormalInit>(1, 0);
    auto x = Param({1}, "x", init);
    auto y = Param({1}, "y", init);
    Minimize(Sqr(x + 2 * y - 7) + Sqr(2 * x + y - 5), 0.01, 1000);

    REQUIRE(abs(x->v[0] - 1) < 1e-5);
    REQUIRE(abs(y->v[0] - 3) < 1e-5);
}

TEST_CASE("diff: minimize rastrigin", "[diff]") {
    auto init = make_shared<NormalInit>(1, 0);
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

TEST_CASE("diff: learn perceptron", "[diff]") {
    // model
    auto x = Data({0, 1}, "x");
    // auto y = Data({0, 1}, "y");
    auto r = Data({0, 1}, "r");

    auto init = make_shared<NormalInit>(1, 0);
    //auto a = Param({1}, "a", init);
    //auto b = Param({1}, "b", init);
    auto c = Param({1}, "c", init);

    auto out = Sigmoid(x + c);
    auto loss = MeanSquareError(out, r);
    loss->name = "loss";
    auto accuracy = Mean(Abs(out - r) < 0.5);
    accuracy->name = "accuracy";

    Print(TopoSort({loss}));

    Model model(loss, accuracy, 4);

    // dataset
    const float A = 1, B = 1, C = 0.4;
    UniformInit gen(-1, 1, 0);
    const int Classes = 2;
    const int SamplesPerClass = 5120;
    const int Samples = Classes * SamplesPerClass;
    vtensor data_x({Samples, 1}, 0);
    // vtensor data_y({Samples, 1}, 0);
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
        // data_y[index] = dy;
        data_r[index] = dr;
        index += 1;
    }
    vector<pair<PDiff, tensor>> dataset = {{x, data_x}, /*{y, data_y},*/ {r, data_r}};

    // train!
    println("train ...");
    Print(model.nodes, true, true);
    for (size_t i = 0; i < 1000; i++) {
        model.Epoch(dataset, 0.01, 0);
        Print(model.nodes, true, true);
    }
}

// Classification:
// learn x + c
// learn a*x + b*y
// learn plane in 2d
// learn hyper plane in Nd
// learn xor (with neurons)
// learn circle in 2d
// learn hyper sphere in Nd

// Regression:
// - multiply two inputs: x * y
// - sqr
// - sqrt
// - sine in 1d
// - sine in 2d : sin(sqrt(x^2 + y^2)

TEST_CASE("diff: learn xor", "[diff]") {}
