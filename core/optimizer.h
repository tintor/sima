#pragma once
#include <core/diff.h>

struct Optimizer {
    float alpha = 0.1;

    virtual void Reset() {}
    virtual void Optimize(span<ParamT*> params);
};

struct Momentum : public Optimizer {
    float momentum = 0.9;

    void Reset() override {
        momentum_exp = 1;
        for (auto& e : ags) EACH(e) e[i] = 0;
    }
    void Optimize(span<ParamT*> params) override;

   private:
    float momentum_exp = 1;
    vector<vtensor> ags;
};

struct RMSProp : public Optimizer {
    float rmsprop = 0.999;
    float rmsprop_eps = 1e-8;

    void Reset() override {
        rmsprop_exp = 1;
        for (auto& e : aggs) EACH(e) e[i] = 0;
    }
    void Optimize(span<ParamT*> params) override;

   private:
    float rmsprop_exp = 1;
    vector<vtensor> aggs;
};

struct Adam : public Optimizer {
    float momentum = 0.9;
    float rmsprop = 0.999;
    float rmsprop_eps = 1e-8;

    void Reset() override {
        rmsprop_exp = 1;
        for (auto& e : ags) EACH(e) e[i] = 0;
        for (auto& e : aggs) EACH(e) e[i] = 0;
    }
    void Optimize(span<ParamT*> params) override;

   private:
    float momentum_exp = 1;
    float rmsprop_exp = 1;
    vector<vtensor> ags, aggs;
};
