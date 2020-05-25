#include <core/diff.h>

thread_local bool Diff::has_overload;

bool IsBroadcastable(tensor_shape a, tensor_shape b) { return a.size != 0 && (a.volume == 1 || a == b.last(a.size)); }

struct BroadcastS : public Diff1 {
    BroadcastS(PDiff a, tensor_shape b) : Diff1(a) { Reshape(b); }
    void Forward() override { EACH(v) v[i] = va[0]; }
    void Backward() override {
        if (ga) EACH(g) ga[0] += g[i];
    }
};

struct BroadcastT : public Diff1 {
    BroadcastT(PDiff a, tensor_shape b) : Diff1(a) { Reshape(b); }
    void Forward() override { EACH(v) v[i] = va[i % va.size]; }
    void Backward() override {
        if (ga) EACH(g) ga[i % ga.size] += g[i];
    }
};

PDiff Broadcast(PDiff a, tensor_shape b) {
    if (a->shape == b) return a;
    if (a->size == 1) return make_shared<BroadcastS>(a, b);
    if (a->shape == b.last(a->rank)) return make_shared<BroadcastT>(a, b);
    Check(false);
    return nullptr;
}

MaxPool2D::MaxPool2D(PDiff a) : Diff1(a) {
    Check(a->rank == 2);
    const uint m = (a->shape(0) + 1) / 2;
    const uint n = (a->shape(1) + 1) / 2;
    Reshape({m, n});
}

void MaxPool2D::Forward() {
    // TODO edge condition on last row / column
    Check(a->shape(0) % 2 == 0);
    Check(a->shape(1) % 2 == 0);

    const uint m = shape(0);
    const uint n = shape(1);
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
    const uint m = shape(0);
    const uint n = shape(1);
    for (uint i = 0; i < m; i++) {
        for (uint j = 0; j < n; j++) {
            const uint p = i * 2, q = j * 2;
            g(i, j) = max(va(p, q), va(p + 1, q), va(p, q + 1), va(p + 1, q + 1));
        }
    }
}
