#pragma once
#include <core/std.h>
#include <core/util.h>

// TODO first dimension can be 0 for mini-batches
struct tensor_shape {
    tensor_shape(uint16_t a = 0, uint16_t b = 0, uint16_t c = 0, uint16_t d = 0) : dim{a, b, c, d} {
        if (d != 0) Check(c != 0);
        if (c != 0) Check(b != 0);
        if (b != 0) Check(a != 0);
    }

    tensor_shape(const tensor_shape& o) : dim{o.dim} {}
    uint16_t operator[](int i) const { return dim[i]; }

    tensor_shape set(int i, uint16_t v) const {
        tensor_shape s = *this;
        s.dim[i] = v;
        return s;
    }

    size_t size() const {
        for (int i = dim.size() - 1; i >= 0; i--) if (dim[i]) return i + 1;
        return 0;
    }

    const uint16_t* begin() const { return dim.data(); }
    const uint16_t* end() const { return dim.data() + size(); }

    size_t volume() const { return Product<size_t>(cspan<uint16_t>(dim.data(), size())); }

    tensor_shape remove_zeros() const {
        tensor_shape s;
        int w = 0;
        for (int i = 0; i < dim.size(); i++) if (dim[i]) s.dim[w++] = dim[i];
        return s;
    }

    tensor_shape pop_front() const {
        tensor_shape s;
        for (int i = 1; i < dim.size(); i++) s.dim[i - 1] = dim[i];
        s.dim[dim.size() - 1] = 0;
        return s;
    }

    tensor_shape pop_back() const {
        tensor_shape s = *this;
        for (int i = dim.size() - 1; i >= 0; i--) if (dim[i]) {
            s.dim[i] = 0;
            break;
        }
        return s;
    }

    auto back() const { auto s = size(); return (s == 0) ? 0 : dim[s - 1]; }

    // first m elements
    tensor_shape first(int m) const {
        tensor_shape s = *this;
        for (int i = m; i < dim.size(); i++) s.dim[i] = 0;
        return s;
    }

    // last m elements
    tensor_shape last(int m) const {
        m = min<int>(m, size());
        tensor_shape s;
        for (int i = 0; i < m; i++) s.dim[i] = dim[i + size() - m];
        for (int i = m; i < dim.size(); i++) s.dim[i] = 0;
        return s;
    }

    tensor_shape push_front(uint16_t a) {
        Check(dim[3] == 0);
        return tensor_shape(a, dim[0], dim[1], dim[2]);
    }

    bool operator==(tensor_shape o) const { return dim == o.dim; }
    bool operator!=(tensor_shape o) const { return dim != o.dim; }

private:
    array<uint16_t, 4> dim;
};

template <typename T = float>
class VTensor;

// Tensor is a pointer to multi-dimensional array (doesn't own memory)
template <typename T = float>
class Tensor {
   public:
    using type = T;

    Tensor() : m_data(nullptr) {}
    Tensor(T* data, tensor_shape shape) : m_data(data), m_shape(shape) {}

    Tensor(const Tensor<T>& o) : m_data(o.m_data), m_shape(o.m_shape) { }

    operator bool() const { return m_data != nullptr; }

    T* data() { return m_data; }
    T const* data() const { return m_data; }

    const T* begin() const { return data(); }
    const T* end() const { return data() + size(); }
    T* begin() { return data(); }
    T* end() { return data() + size(); }

    size_t size() const { return m_shape.volume(); }
    T& operator[](size_t index) { return m_data[index]; }
    const T& operator[](size_t index) const { return m_data[index]; }

    const auto& shape() const { return m_shape; }
    T& operator()(uint16_t a, uint16_t b) { return operator()({a, b}); }
    T& operator()(uint16_t a, uint16_t b, uint16_t c) { return operator()({a, b, c}); }
    T& operator()(uint16_t a, uint16_t b, uint16_t c, uint16_t d) { return operator()({a, b, c, d}); }
    T& operator()(tensor_shape index) { return m_data[offset(index)]; }
    const T& operator()(tensor_shape index) const { return m_data[offset(index)]; }

    size_t offset(tensor_shape index) const {
        size_t out = index[0];
        for (size_t i = 1; i < index.size(); i++) out = out * m_shape[i] + index[i];
        return out;
    }

    // sub-tensor for the fixed value of the first dimension
    Tensor slice(uint16_t a) {
        auto v = m_shape.pop_front().volume();
        return Tensor(m_data + a * v, v);
    }

    const Tensor slice(uint16_t a) const {
        auto v = m_shape.pop_front().volume();
        return Tensor(m_data + a * v, v);
    }

   private:
    T* m_data;
    tensor_shape m_shape;
};

// VTensor is multi-dimensional vector (owns memory)
// Can't be reshaped.
template <typename T>
class VTensor {
   public:
    VTensor() {}

    VTensor(tensor_shape shape, T init = T()) : m_data(Product<size_t>(shape), init), m_shape(shape) {}

    VTensor(const Tensor<T>& o) { operator=(o); }

    VTensor(const VTensor<T>& o) { operator=(o); }

    operator bool() const { return m_data.size() > 0; }

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

    T& operator()(uint16_t a, uint16_t b) { return operator()({a, b}); }
    T& operator()(uint16_t a, uint16_t b, uint16_t c) { return operator()({a, b, c}); }
    T& operator()(uint16_t a, uint16_t b, uint16_t c, uint16_t d) { return operator()({a, b, c, d}); }
    T& operator()(tensor_shape index) { return m_data[offset(index)]; }

    const T& operator()(uint16_t a, uint16_t b) const { return operator()({a, b}); }
    const T& operator()(uint16_t a, uint16_t b, uint16_t c) const { return operator()({a, b, c}); }
    const T& operator()(uint16_t a, uint16_t b, uint16_t c, uint16_t d) const { return operator()({a, b, c, d}); }
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
