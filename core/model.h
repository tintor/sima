#pragma once
#include <core/optimizer.h>

// Returns a list of all diffs that need to be computed (in order) before all goal diffs.
// Goal diffs are also ordered if one depends on any other.
vector<PDiff> TopoSort(const cspan<PDiff> heads);

// TODO
inline vector<PDiff> LoadModel(string_view filename) { return {}; }

// TODO - save model structure and v tensor of all Param nodes
inline void SaveModel(span<PDiff> model, string_view filename) {}

using Metrics = unordered_map<string, float>;

class Model {
   public:
    shared_ptr<Optimizer> optimizer = make_shared<Optimizer>();

    Model(PDiff out, PDiff loss, PDiff accuracy);

    void Forward();
    void Backward();
    void Print() const;

    void Iterate(size_t iterations) {
        for (size_t i : range(iterations)) {
            Forward();
            Backward();
        }
    }

    Metrics Epoch(cspan<pair<PDiff, tensor>> data, std::mt19937_64& random, bool verbose = true, uint epoch = 0);

    tensor::type Loss() const { return m_loss->v[0]; }
    tensor::type Accuracy() const { return m_accuracy->v[0]; }

   private:
    const PDiff m_out;
    const PDiff m_loss;
    const PDiff m_accuracy;
    vector<PDiff> m_nodes;
    vector<Diff*> m_forward_nodes, m_backward_nodes, m_end_epoch_nodes;
    vector<ParamT*> m_params;
    vector<uint> m_samples;
};

inline void Minimize(PDiff loss, float alpha, size_t iterations) {
    Model model(nullptr, loss, nullptr);
    model.optimizer->alpha = alpha;
    model.Iterate(iterations);
}

// subtract mean and divide by stdev
void NormalizeDataset(tensor);
