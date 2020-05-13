#pragma once
#include <core/std.h>
#include <core/util.h>

// TODO:
// Optimized Tensor and TensorView class layout (main difference is that TensorView doesn't own its memory):
// - m_data ptr (based on raw_array)
// - 8 byte shape array: std::array<uint16_t, 4> (if >4 dimensions are needed then shape vector can be varint encoded)
// simple!
// efficient for small tensors
// supports large tensors
// can be optimized for inline storage for tensor of size 1 float (ie. scalars can be efficient as tensors)

template <typename T>
class Tensor;

template <typename T>
class TensorSpan {
   public:
    TensorSpan(T* data, cspan<uint32_t> shape) : m_data(data), m_shape(shape.size()) { Copy(shape, m_shape); }

    template <typename X>
    TensorSpan(const TensorSpan<X>& o) : m_data(o.m_data) {
        Copy(o.m_shape, m_shape);
    }
    /*template <typename X>
    TensorSpan(const Tensor<X>& o) : m_data(o.data()) {
        Copy(o.shape(), m_shape);
    }*/

    T* data() { return m_data; }
    T const* data() const { return m_data; }

    const T* begin() const { return data(); }
    const T* end() const { return data() + size(); }
    T* begin() { return data(); }
    T* end() { return data() + size(); }

    size_t size() const { return Product<size_t>(m_shape); }
    T& operator[](size_t index) { return m_data[index]; }
    const T& operator[](size_t index) const { return m_data[index]; }

    const auto& shape() const { return m_shape; }
    T& operator()(cspan<uint32_t> index) { return m_data[offset(index)]; }
    const T& operator()(cspan<uint32_t> index) const { return m_data[offset(index)]; }

    size_t offset(cspan<uint32_t> index) const {
        size_t out = index[0];
        for (size_t i = 1; i < index.size(); i++) out = out * m_shape[i] + index[i];
        return out;
    }

   private:
    T* m_data;
    vector<uint32_t> m_shape;
};

template <typename T>
class Tensor {
   public:
    Tensor(cspan<uint32_t> shape, T init = T()) : m_data(Product<size_t>(shape), init), m_shape(shape.size()) {
        Copy(shape, m_shape);
    }

    template <typename X>
    Tensor(const TensorSpan<X>& o) {
        operator=(o);
    }

    Tensor(const Tensor<const T>& o) { operator=(o); }

    T* data() { return m_data.data(); }
    T const* data() const { return m_data.data(); }

    const T* begin() const { return data(); }
    const T* end() const { return data() + size(); }
    T* begin() { return data(); }
    T* end() { return data() + size(); }

    size_t size() const { return m_data.size(); }
    T& operator[](size_t index) { return m_data[index]; }
    const T& operator[](size_t index) const { return m_data[index]; }

    cspan<uint32_t> shape() const { return m_shape; }
    T& operator()(cspan<uint32_t> index) { return m_data[offset(index)]; }
    const T& operator()(cspan<uint32_t> index) const { return m_data[offset(index)]; }

    operator TensorSpan<T>() { return TensorSpan<T>(m_data.data(), m_shape); }
    operator TensorSpan<const T>() const { return TensorSpan<const T>(m_data.data(), m_shape); }

    TensorSpan<T> sub(uint32_t index) {
        cspan<uint32_t> sub_shape = m_shape.pop_back();
        return TensorSpan<T>(m_data + index * Product<size_t>(sub_shape), sub_shape);
    }
    TensorSpan<const T> sub(uint32_t index) const {
        cspan<uint32_t> sub_shape = cspan<uint32_t>(m_shape).pop_back();
        return TensorSpan<const T>(m_data.data() + index * Product<size_t>(sub_shape), sub_shape);
    }

    size_t offset(cspan<uint32_t> index) const {
        size_t out = index[0];
        for (size_t i = 1; i < index.size(); i++) out = out * m_shape[i] + index[i];
        return out;
    }

    void operator=(const TensorSpan<const float>& o) {
        m_shape.resize(o.shape().size());
        Copy(o.shape(), m_shape);
        m_data.resize(o.size());
        // TODO(Marko) not safe when strides are added to TensorSpan
        std::copy(o.data(), o.data() + o.size(), m_data.data());
    }

    void reshape(cspan<uint32_t> shape, T init = T()) {
        m_shape.resize(shape.size());
        Copy(shape, m_shape);
        m_data.resize(Product<size_t>(m_shape), init);
    }

   private:
    vector<T> m_data;
    vector<uint32_t> m_shape;
};
