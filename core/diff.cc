#include <core/diff.h>

thread_local bool Diff::has_overload;

bool IsBroadcastable(dim4 a, dim4 b) { return a.ndims() != 0 && (a.elements() == 1 /*|| a == b.last(a.size)*/); }

struct BroadcastS : public Diff1 {
    BroadcastS(PDiff a, dim4 b) : Diff1(a) { Reshape(b); }
    void Forward(bool) override { EACH(v) v[i] = va[0]; }
    void Backward() override {
        if (ga) EACH(g) ga[0] += g[i];
    }
};

/*struct BroadcastT : public Diff1 {
    BroadcastT(PDiff a, dim4 b) : Diff1(a) { Reshape(b); }
    void Forward() override { EACH(v) v[i] = va[i % va.size]; }
    void Backward() override {
        if (ga) EACH(g) ga[i % ga.size] += g[i];
    }
};*/

PDiff Broadcast(PDiff a, dim4 b) {
    if (a->shape() == b) return a;
    if (a->elements() == 1) return make_shared<BroadcastS>(a, b);
    // if (a->shape == b.last(a->rank)) return make_shared<BroadcastT>(a, b);
    Check(false);
    return nullptr;
}

MaxPool2D::MaxPool2D(PDiff a) : Diff1(a) {
    Check(a->ndims() == 2);
    const uint m = (a->dim(0) + 1) / 2;
    const uint n = (a->dim(1) + 1) / 2;
    Reshape({m, n});
}

void MaxPool2D::Forward(bool) {
    // TODO edge condition on last row / column
    Check(a->dim(0) % 2 == 0);
    Check(a->dim(1) % 2 == 0);

    const uint m = dim(0);
    const uint n = dim(1);
    for (uint i = 0; i < m; i++) {
        for (uint j = 0; j < n; j++) {
            const uint p = i * 2, q = j * 2;
            v(i, j) = max(va(p, q), va(p + 1, q), va(p, q + 1), va(p + 1, q + 1));
        }
    }
}

void MaxPool2D::Backward() {
    if (!ga) return;
    Check(false);
    const uint m = dim(0);
    const uint n = dim(1);
    for (uint i = 0; i < m; i++) {
        for (uint j = 0; j < n; j++) {
            const uint p = i * 2, q = j * 2;
            g(i, j) = max(va(p, q), va(p + 1, q), va(p, q + 1), va(p + 1, q + 1));
        }
    }
}
