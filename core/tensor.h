#pragma once
#include <core/format.h>
#include <core/std.h>
#include <core/util.h>

struct tensor_shape {
    tensor_shape(uint16_t a = 0, uint16_t b = 0, uint16_t c = 0, uint16_t d = 0) : dim{a, b, c, d} {}

    tensor_shape(const tensor_shape& o) : dim{o.dim} {}
    uint16_t operator[](int i) const { return dim[i]; }

    tensor_shape set(int i, uint16_t v) const {
        tensor_shape s = *this;
        s.dim[i] = v;
        return s;
    }

    size_t size() const {
        for (int i = dim.size() - 1; i >= 0; i--)
            if (dim[i]) return i + 1;
        return 0;
    }

    size_t volume() const { return (size() == 0) ? 0 : Product<size_t>(cspan<uint16_t>(dim.data(), size())); }

    tensor_shape remove_zeros() const {
        tensor_shape s;
        int w = 0;
        for (int i = 0; i < dim.size(); i++)
            if (dim[i]) s.dim[w++] = dim[i];
        return s;
    }

    tensor_shape remove_ones() const {
        tensor_shape s;
        int w = 0;
        for (int i = 0; i < dim.size(); i++)
            if (dim[i] != 1) s.dim[w++] = dim[i];
        return s;
    }

    tensor_shape concat(tensor_shape o) const {
        Check(size() + o.size() <= dim.size());
        tensor_shape s = *this;
        auto m = size();
        for (int i = 0; i < o.size(); i++) s.dim[m++] = o[i];
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
        for (int i = dim.size() - 1; i >= 0; i--)
            if (dim[i]) {
                s.dim[i] = 0;
                break;
            }
        return s;
    }

    auto back() const {
        auto s = size();
        return (s == 0) ? 0 : dim[s - 1];
    }

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

    tensor_shape push_front(uint16_t a) const {
        Check(dim[3] == 0);
        return tensor_shape(a, dim[0], dim[1], dim[2]);
    }

    tensor_shape push_back(uint16_t a) const {
        Check(dim[3] == 0);
        tensor_shape s = *this;
        s.dim[s.size()] = a;
        return s;;
    }

    operator string() const {
        ostringstream os;
        os << '[';
        for (int i = 0; i < size(); i++) {
            if (i > 0) os << ' ';
            os << dim[i];
        }
        os << ']';
        return os.str();
    }

    bool operator==(tensor_shape o) const { return dim == o.dim; }
    bool operator!=(tensor_shape o) const { return dim != o.dim; }

   private:
    array<uint16_t, 4> dim;
};

// Tensor is a pointer to multi-dimensional array (doesn't own memory)
template <typename T = float>
class Tensor {
   public:
    using type = T;

    Tensor() : m_data(nullptr) {}
    Tensor(T* data, tensor_shape shape) : m_data(data), m_shape(shape) {
        Check(data || shape.volume() == 0);
        Check(shape.volume() != 0 || data == nullptr);
    }

    Tensor(const Tensor<T>& o) : m_data(o.m_data), m_shape(o.m_shape) {}

    operator bool() const { return m_data != nullptr; }

    T* data() { return m_data; }
    T const* data() const { return m_data; }

    const T* begin() const { return data(); }
    const T* end() const { return data() + size(); }
    T* begin() { return data(); }
    T* end() { return data() + size(); }

    size_t size() const { return m_shape.volume(); }
    T& operator[](size_t index) {
        DCheck(index < size(), format("%s < %s", index, size()));
        return m_data[index];
    }
    const T& operator[](size_t index) const {
        DCheck(index < size(), format("%s < %s", index, size()));
        return m_data[index];
    }

    const auto& shape() const { return m_shape; }
    T& operator()(uint16_t a, uint16_t b) { return operator()({a, b}); }
    T& operator()(uint16_t a, uint16_t b, uint16_t c) { return operator()({a, b, c}); }
    T& operator()(uint16_t a, uint16_t b, uint16_t c, uint16_t d) { return operator()({a, b, c, d}); }
    T& operator()(tensor_shape index) { return m_data[offset(index)]; }
    const T& operator()(tensor_shape index) const { return operator[](offset(index)); }

    size_t offset(tensor_shape index) const {
        size_t out = index[0];
        for (size_t i = 1; i < m_shape.size(); i++) out = out * m_shape[i] + index[i];
        return out;
    }

    // sub-tensor for the fixed value of the first dimension
    Tensor slice(uint16_t a) {
        DCheck(m_shape.size() > 0 && a < m_shape[0], "");
        auto v = m_shape.pop_front().volume();
        return Tensor(m_data + a * v, v);
    }

    const Tensor slice(uint16_t a) const {
        DCheck(m_shape.size() > 0 && a < m_shape[0], "");
        auto v = m_shape.pop_front().volume();
        return Tensor(m_data + a * v, v);
    }

    void copy_from(const Tensor& o) {
        Check(shape() == o.shape());
        std::copy(o.begin(), o.end(), begin());
    }

    bool operator==(const Tensor& o) const {
        return m_shape == o.m_shape && std::equal(m_data, m_data + size(), o.m_data);
    }
    bool operator!=(const Tensor& o) const { return !operator==(o); }

    Tensor& operator=(const Tensor& o) {
        m_data = o.m_data;
        m_shape = o.m_shape;
        return *this;
    }

    operator string() const {
        string s;
        s += '[';
        for (size_t i = 0; i < size(); i++) {
            if (i > 0) s += ' ';
            format_s(s, "%s", m_data[i]);
        }
        s += ']';
        return s;
    }

   protected:
    T* m_data;
    tensor_shape m_shape;
};

// VTensor is multi-dimensional vector (owns memory)
template <typename T>
class VTensor : public Tensor<T> {
   public:
    VTensor() {}
    VTensor(tensor_shape shape, T init = T()) : m_vector(shape.volume(), init) {
        Tensor<T>::m_shape = shape;
        Tensor<T>::m_data = m_vector.data();
    }
    VTensor(VTensor<T>&& o) { operator=(o); }
    VTensor(const Tensor<T>& o) { operator=(o); }

    VTensor& operator=(VTensor<T>&& o) {
        if (this == &o) return *this;
        Tensor<T>::m_shape = o.shape();
        m_vector.clear();
        swap(m_vector, o.m_vector);
        Tensor<T>::m_data = m_vector.data();

        o.m_data = nullptr;
        o.m_shape = tensor_shape();
        return *this;
    }

    VTensor& operator=(const Tensor<T> o) {
        if (this == &o) return *this;
        Tensor<T>::m_shape = o.shape();
        m_vector.resize(o.size());
        Tensor<T>::m_data = m_vector.data();
        // TODO not safe if strides are added to Tensor
        std::copy(o.data(), o.data() + o.size(), m_vector.data());
        return *this;
    }

    void reshape(tensor_shape new_shape, T init = T()) {
        Tensor<T>::m_shape = new_shape;
        m_vector.resize(new_shape.volume(), init);
        Tensor<T>::m_data = m_vector.data();
    }

   private:
    vector<T> m_vector;
};

using tensor = Tensor<float>;
using vtensor = VTensor<float>;
