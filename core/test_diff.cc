#include <core/model.h>
#include <core/thread.h>

#include <catch.hpp>

TEST_CASE("diff: minimize circle", "[diff]") {
    auto x = Param({1}) << "x";
    auto y = Param({1}) << "y";
    x->v[0] = 3.14f;
    y->v[0] = 2.16f;
    Minimize(Sqr(x) + Sqr(y), 0.1, 70);

    REQUIRE(abs(x->v[0]) < 1e-6);
    REQUIRE(abs(y->v[0]) < 1e-6);
}

void Test(string_view name, Model& model, PDiff x, PDiff ref, PDiff b, vector<string>& table) {
    auto row = table.begin();
    if (name == "sgd") format_s(*row, "b|");
    format_s(*row++, "%s|", name);

    for (float e : range(-1.f, 1.01f, 0.1f)) {
        x->v[0] = e;
        ref->v[0] = 0;
        model.optimizer->Reset();
        model.Forward();
        if (name == "sgd") format_s(*row, "%s|", b->v[0]);
        int count = 250;
        for (int i = 0; i < 250; i++) {
            model.Forward();
            model.Backward();
            if (b->v[0] <= 0.05) {
                count = i + 1;
                break;
            }
        }
        REQUIRE(b->v[0] <= 0.05);
        format_s(*row++, "%s|", count);
    }

    for (float e : range(-1.f, 1.01f, 0.1f)) {
        x->v[0] = e;
        ref->v[0] = 1;
        model.optimizer->Reset();
        model.Forward();
        if (name == "sgd") format_s(*row, "%s|", b->v[0]);
        int count = 250;
        for (int i = 0; i < 250; i++) {
            model.Forward();
            model.Backward();
            if (b->v[0] >= 0.95) {
                count = i + 1;
                break;
            }
        }
        REQUIRE(b->v[0] >= 0.95);
        format_s(*row++, "%s|", count);
    }
}

TEST_CASE("diff: minimize binary cross entropy", "[diff]") {
    auto x = Param({1}) << "x";
    auto ref = Const(0);
    auto b = Logistic(x) << "b";
    auto loss = BinaryCrossEntropy<true>(ref, b);
    Model model(loss);

    vector<string> table;
    table << format("x|r|");
    for (float e : range(-1.f, 1.01f, 0.1f)) table << format("%+.1f|0|", e);
    for (float e : range(-1.f, 1.01f, 0.1f)) table << format("%+.1f|1|", e);

    Test("sgd", model, x, ref, b, table);

    model.optimizer = make_shared<Momentum>();
    Test("momentum", model, x, ref, b, table);

    model.optimizer = make_shared<RMSProp>();
    Test("rmsprop", model, x, ref, b, table);

    model.optimizer = make_shared<Adam>();
    Test("adam", model, x, ref, b, table);

    PrintTable(table, '|', " ", {3, 4, 5, 6});
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

TEST_CASE("diff: learn perceptron, plane in 2d", "[diff]") {
    constexpr bool D3 = false;
    parallel(5, [&](size_t seed) {
        const int Batch = 24;
        auto x = Data({Batch, 1}) << "x";
        auto y = Data({Batch, 1}) << "y";
        auto z = Data({Batch, 1}) << "y";
        auto ref = Data({Batch, 1}) << "ref";

        std::mt19937_64 random(seed);
        auto init = make_shared<NormalInit>(1 / sqrt(2), random);
        auto a = Param({1}, init) << "a";
        auto b = Param({1}, init) << "b";
        auto c = Param({1}) << "d";
        auto d = D3 ? Param({1}, init) << "c" : nullptr;

        auto fc = x * a + y * b + c;
        if (D3) fc = fc + z * d;
        auto out = Logistic(fc, 15) << "out";
        auto loss = BinaryCrossEntropy(ref, out) << "loss";

        Model model(loss, Accuracy(ref, out));
        model.optimizer = make_shared<Adam>();
        model.optimizer->alpha = 0.1;

        // dataset
        const float A = 0.2, B = 0.4, C = -0.8, D = 0.1;
        UniformInit gen(-1, 1, random);
        const int Classes = D3 ? 3 : 2;
        const int SamplesPerClass = 20000 - 8;
        const int Samples = Classes * SamplesPerClass;
        vtensor data_x({Samples, 1}, 0);
        vtensor data_y({Samples, 1}, 0);
        vtensor data_z({Samples, 1}, 0);
        vtensor data_r({Samples, 1}, 0);
        int count[2] = {0, 0};
        int index = 0;
        while (index < Samples) {
            float dx = gen.get();
            float dy = gen.get();
            float dz = D3 ? gen.get() : 0;
            float dr = (dx * A + dy * B + dz * C + D >= 0) ? 1 : 0;
            int c = round(dr);
            if (count[c] >= SamplesPerClass) continue;
            count[c] += 1;
            data_x[index] = dx;
            data_y[index] = dy;
            if (D3) data_z[index] = dz;
            data_r[index] = dr;
            index += 1;
        }
        vector<pair<PDiff, tensor>> dataset;
        if (D3) dataset = {{x, data_x}, {y, data_y}, {z, data_z}, {ref, data_r}};
        if (!D3) dataset = {{x, data_x}, {y, data_y}, {ref, data_r}};
        NormalizeDataset(data_x);
        NormalizeDataset(data_y);
        if (D3) NormalizeDataset(data_z);

        // train!
        Metrics metrics;
        for (auto i : range(1000)) metrics = model.Epoch(dataset, random, false, i);
        REQUIRE(metrics.at("accuracy") >= 0.9994);
        println("accuracy: %s", metrics.at("accuracy"));
    });
}

PDiff Neuron(PDiff x, PDiff y, shared_ptr<Init> init) {
    auto a = Param({1}, init) << "a";
    auto b = Param({1}, init) << "b";
    auto c = Param({1}) << "c";
    return Logistic(x * a + y * b + c, 15);
}

PDiff Neuron(PDiff x, PDiff y, PDiff z, shared_ptr<Init> init) {
    auto a = Param({1}, init) << "a";
    auto b = Param({1}, init) << "b";
    auto c = Param({1}, init) << "c";
    auto d = Param({1}) << "d";
    return Logistic(x * a + y * b + z * c + d, 15);
}

// Effect of batch size on duration and accuracy
// B16 -> 11.7s 99.17%
// B20 -> 10.9s 99.23%
// B24 -> 10.5s 99.23%
// B32 -> 8.7s 99.19%
// B64 -> 7.8s 98.98%
// B128 -> 7.2s 98.43%
// B256 -> 7.0s 97.39%
// B512 -> 6.9s 95.28%
TEST_CASE("diff: learn two layer network, circle in 2d", "[diff]") {
    parallel(5, [&](size_t seed) {
        const int Batch = 24;
        auto x = Data({Batch, 1}) << "x";
        auto y = Data({Batch, 1}) << "y";
        auto ref = Data({Batch, 1}) << "ref";

        std::mt19937_64 random(seed);
        auto init = make_shared<NormalInit>(1, random);
        auto h1 = Neuron(x, y, init) << "h1";
        auto h2 = Neuron(x, y, init) << "h2";
        auto h3 = Neuron(x, y, init) << "h3";
        auto init2 = make_shared<NormalInit>(1, random);
        auto out = Neuron(h1, h2, h3, init2) << "out";
        auto loss = MeanSquareError(ref, out) << "loss";
        Model model(loss, Accuracy(ref, out));
        model.optimizer->alpha = 0.1;

        // dataset
        UniformInit gen(-1, 1, random);
        const int Classes = 2;
        const int SamplesPerClass = 20000 - 8;
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
        Metrics metrics;
        for (auto i : range(1000)) metrics = model.Epoch(dataset, random, false, i);
        println("accuracy: %s", metrics.at("accuracy"));
        REQUIRE(metrics.at("accuracy") >= 0.9922);
    });
}

TEST_CASE("diff: learn FC perceptron, hyperplane", "[diff]") {
    constexpr int N = 2;
    parallel(5, [&](size_t seed) {
        const int Batch = 24;
        auto in = Data({Batch, N}) << "in";
        auto ref = Data({Batch, 1}) << "ref";

        std::mt19937_64 random(seed);
        auto init = make_shared<NormalInit>(1, random);
        auto fc = FullyConnected(in, 1, init);
        auto out = Logistic(fc, 15) << "out";
        auto loss = BinaryCrossEntropy(ref, out) << "loss";

        Model model(loss, Accuracy(ref, out));
        model.optimizer = make_shared<Adam>();
        model.optimizer->alpha = 0.1;

        // dataset
        UniformInit gen(-1, 1, random);
        vector<float> W;
        for (int i : range(N)) W.push_back(gen.get() / 2);
        float B = gen.get() / 10;

        const int Classes = 2;
        const int SamplesPerClass = 20000 - 8;
        const int Samples = Classes * SamplesPerClass;
        vtensor data_in({Samples, N}, 0);
        vtensor data_ref({Samples, 1}, 0);
        int count[2] = {0, 0};
        int index = 0;
        vector<float> d;
        while (index < Samples) {
            d.clear();
            for (int i : range(N)) d.push_back(gen.get());
            float s = B;
            for (int i : range(N)) s += d[i] * W[i];
            float dr = (s >= 0) ? 1 : 0;
            int c = round(dr);
            if (count[c] >= SamplesPerClass) continue;
            count[c] += 1;
            for (int i : range(N)) data_in(index, i) = d[i];
            data_ref[index] = dr;
            index += 1;
        }
        vector<pair<PDiff, tensor>> dataset = {{in, data_in}, {ref, data_ref}};
        NormalizeDataset(data_in);

        // train!
        Metrics metrics;
        for (auto i : range(1000)) metrics = model.Epoch(dataset, random, false, i);
        println("accuracy: %s", metrics.at("accuracy"));
        REQUIRE(metrics.at("accuracy") >= 0.09994);
    });
}

TEST_CASE("diff: learn FC two layer network, circle in 2d", "[diff]") {
    parallel(5, [&](size_t seed) {
        const int Batch = 24;
        auto in = Data({Batch, 2}) << "in";
        auto ref = Data({Batch, 1}) << "ref";

        std::mt19937_64 random(seed);
        auto init = make_shared<NormalInit>(1, random);

        auto x = in;
        x = FullyConnected(x, 3, init) << "fc";
        x = Logistic(x, 15);
        x = FullyConnected(x, 1, init);
        x = Logistic(x, 15);

        auto out = x << "out";
        auto loss = MeanSquareError(ref, out) << "loss";

        Model model(loss, Accuracy(ref, out));
        model.optimizer->alpha = 0.1;

        // dataset
        UniformInit gen(-1, 1, random);
        const int Classes = 2;
        const int SamplesPerClass = 20000 - 8;
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
        NormalizeDataset(data_in);

        // train!
        Metrics metrics;
        for (auto i : range(1000)) metrics = model.Epoch(dataset, random, false, i);
        println("accuracy: %s", metrics.at("accuracy"));
        REQUIRE(metrics.at("accuracy") >= 0.9617);
    });
}

// Classification:
// learn hyper plane in Nd
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
