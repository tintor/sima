#pragma once
#include <core/std.h>
#include <core/util.h>

// TODO:
// Optimized VTensor and VTensorView class layout (main difference is that VTensorView doesn't own its memory):
// - m_data ptr (based on raw_array)
// - 8 byte shape array: std::array<uint16_t, 4> (if >4 dimensions are needed then shape vector can be varint encoded)
// simple!
// efficient for small Tensors
// supports large Tensors
// can be optimized for inline storage for Tensor of size 1 float (ie. scalars can be efficient as Tensors)

using tensor_shape = array<uint16_t, 4>;

inline tensor_shape Shape(uint16_t a) {
    return {a, 0, 0, 0};
}

inline tensor_shape Shape(uint16_t a, uint16_t b) {
    return {a, b, 0, 0};
}

inline tensor_shape Concat(uint16_t a, tensor_shape b) {
    Check(b[3] == 0);
    return {a, b[0], b[1], b[2]};
}

inline tensor_shape PopBack(tensor_shape a) {
    for (int i = 3; i >= 0; i--)
        if (a[i] != 0) {
            a[i] = 0;
            return a;
        }
    Check(false);
    return {};
}

template <typename T = float>
class VTensor;

// Tensor is a pointer to multi-dimensional array (doesn't own memory)
template <typename T = float>
class Tensor {
   public:
    using type = T;

    Tensor(T* data, tensor_shape shape) : m_data(data), m_shape(shape) {}

    Tensor(const Tensor<T>& o) : m_data(o.m_data) { Copy(o.m_shape, m_shape); }

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
    T& operator()(tensor_shape index) { return m_data[offset(index)]; }
    const T& operator()(tensor_shape index) const { return m_data[offset(index)]; }

    size_t offset(tensor_shape index) const {
        size_t out = index[0];
        for (size_t i = 1; i < index.size(); i++) out = out * m_shape[i] + index[i];
        return out;
    }

   private:
    T* m_data;
    tensor_shape m_shape;
};

// VTensor is multi-dimensional vector (owns memory)
template <typename T>
class VTensor {
   public:
    VTensor(tensor_shape shape, T init = T()) : m_data(Product<size_t>(shape), init), m_shape(shape) {}

    VTensor(const Tensor<T>& o) { operator=(o); }

    VTensor(const VTensor<T>& o) { operator=(o); }

    T* data() { return m_data.data(); }
    T const* data() const { return m_data.data(); }

    const T* begin() const { return data(); }
    const T* end() const { return data() + size(); }
    T* begin() { return data(); }
    T* end() { return data() + size(); }

    size_t size() const { return m_data.size(); }
    T& operator[](size_t index) { return m_data[index]; }
    const T& operator[](size_t index) const { return m_data[index]; }

    tensor_shape shape() const { return m_shape; }
    T& operator()(tensor_shape index) { return m_data[offset(index)]; }
    const T& operator()(tensor_shape index) const { return m_data[offset(index)]; }

    operator Tensor<T>() { return Tensor<T>(m_data.data(), m_shape); }
    operator const Tensor<T>() const { return Tensor<T>(const_cast<T*>(m_data.data()), m_shape); }

#if 0
    Tensor<T> sub(uint32_t index) {
        tensor_shape sub_shape = m_shape.pop_back();
        return Tensor<T>(m_data + index * Product<size_t>(sub_shape), sub_shape);
    }
    Tensor<const T> sub(uint32_t index) const {
        tensor_shape sub_shape = tensor_shape(m_shape).pop_back();
        return Tensor<const T>(m_data.data() + index * Product<size_t>(sub_shape), sub_shape);
    }
#endif

    size_t offset(tensor_shape index) const {
        size_t out = index[0];
        for (size_t i = 1; i < index.size(); i++) out = out * m_shape[i] + index[i];
        return out;
    }

    void operator=(const Tensor<float> o) {
        m_shape = o.shape();
        m_data.resize(o.size());
        // TODO(Marko) not safe when strides are added to Tensor
        std::copy(o.data(), o.data() + o.size(), m_data.data());
    }

    void reshape(tensor_shape shape, T init = T()) {
        m_shape = shape;
        m_data.resize(Product<size_t>(m_shape), init);
    }

   private:
    vector<T> m_data;
    tensor_shape m_shape;
};

using tensor = Tensor<float>;
using vtensor = VTensor<float>;
