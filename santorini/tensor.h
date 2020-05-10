#pragma once
#include <core/std.h>

template <typename T, typename Range>
T Product(const Range& range) {
    T out = 1;
    for (const auto& e : range) out *= e;
    return out;
}

template <typename CollectionA, typename CollectionB>
void Copy(const CollectionA& a, CollectionB& b) {
    b.resize(a.size());
    for (size_t i = 0; i < a.size(); i++) b[i] = a[i];
}

template <typename T>
class Tensor {
   public:
    Tensor(cspan<uint32_t> shape) : m_shape(shape.size()), m_data(Product<uint32_t>(shape)) { Copy(shape, m_shape); }

    size_t size() const { return m_data.size(); }
    const auto& shape() const { return m_shape; }

    T& operator[](size_t index) { return m_data[index]; }
    const T& operator[](size_t index) const { return m_data[index]; }

    T& operator()(cspan<uint32_t> index) { return m_data[offset(index)]; }
    const T& operator()(cspan<uint32_t> index) const { return m_data[offset(index)]; }

    size_t offset(cspan<uint32_t> index) const {
        size_t out = index[0];
        for (size_t i = 1; i < index.size(); i++) out = out * m_shape[i] + index[i];
        return out;
    }

   private:
    vector<T> m_data;
    vector<uint32_t> m_shape;
};
