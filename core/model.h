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

    Model() {}
    explicit Model(std::initializer_list<PDiff> heads);

    void Forward();
    void Backward(PDiff loss);
    void Print() const;

    void Iterate(size_t iterations, PDiff loss) {
        for (size_t i : range(iterations)) {
            Forward();
            Backward(loss);
        }
    }

    Metrics Epoch(PDiff loss, PDiff accuracy, cspan<pair<PDiff, tensor>> data, std::mt19937_64& random, bool verbose = true, uint epoch = 0);

   private:
    vector<PDiff> m_nodes;
    vector<Diff*> m_forward_nodes, m_backward_nodes, m_end_epoch_nodes;
    vector<ParamT*> m_params;
    vector<uint> m_samples;
};

inline void Minimize(PDiff loss, float alpha, size_t iterations) {
    Model model({loss});
    model.optimizer->alpha = alpha;
    model.Iterate(iterations, loss);
}

// subtract mean and divide by stdev
void NormalizeDataset(tensor);
