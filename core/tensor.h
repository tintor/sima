#pragma once
#include <core/format.h>
#include <core/std.h>
#include <core/util.h>

struct tensor_shape {
    tensor_shape() {
        for (size_t i = 0; i < dim.size(); i++)
            dim[i] = 0;
    }

    tensor_shape(std::initializer_list<uint> d) {
        for (size_t i = 0; i < dim.size(); i++)
            dim[i] = (i < d.size()) ? d.begin()[i] : 0;
    }

    tensor_shape(const tensor_shape& o) : dim{o.dim} {}
    uint operator[](int i) const { return dim[i]; }

    tensor_shape set(int i, uint v) const {
        tensor_shape s = *this;
        s.dim[i] = v;
        return s;
    }

    size_t size() const {
        for (int i = dim.size() - 1; i >= 0; i--)
            if (dim[i]) return i + 1;
        return 0;
    }

    size_t volume() const { return (size() == 0) ? 0 : Product<size_t>(cspan<uint>(dim.data(), size())); }

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

    uint back() const {
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

    tensor_shape push_front(uint a) const {
        Check(dim[3] == 0);
        return tensor_shape{a, dim[0], dim[1], dim[2]};
    }

    tensor_shape push_back(uint a) const {
        Check(dim[3] == 0);
        tensor_shape s = *this;
        s.dim[s.size()] = a;
        return s;
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
    array<uint, 4> dim;
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

    T* begin() { return data(); }
    T* end() { return data() + size(); }
    const T* begin() const { return data(); }
    const T* end() const { return data() + size(); }

    size_t size() const { return m_shape.volume(); }
    const auto& shape() const { return m_shape; }

    T& operator[](size_t index) {
        DCheck(index < size(), format("%s < %s", index, size()));
        return m_data[index];
    }
    const T& operator[](size_t index) const {
        DCheck(index < size(), format("%s < %s", index, size()));
        return m_data[index];
    }

    T& operator()(size_t a) { return m_data[offset(a)]; }
    const T& operator()(size_t a) const { return m_data[offset(a)]; }

    T& operator()(size_t a, size_t b) { return m_data[offset(a, b)]; }
    const T& operator()(size_t a, size_t b) const { return m_data[offset(a, b)]; }

    T& operator()(size_t a, size_t b, size_t c) { return m_data[offset(a, b, c)]; }
    const T& operator()(size_t a, size_t b, size_t c) const { return m_data[offset(a, b, c)]; }

    T& operator()(size_t a, size_t b, size_t c, size_t d) { return m_data[offset(a, b, c, d)]; }
    const T& operator()(size_t a, size_t b, size_t c, size_t d) const { return m_data[offset(a, b, c, d)]; }

    size_t offset(size_t a) const {
        DCheck(m_shape.size() == 1, "");
        DCheck(a < m_shape[0], "");
        return a;
    }
    size_t offset(size_t a, size_t b) const {
        DCheck(m_shape.size() == 2, "");
        DCheck(a < m_shape[0], "");
        DCheck(b < m_shape[1], "");
        return a * m_shape[1] + b;
    }
    size_t offset(size_t a, size_t b, size_t c) const {
        DCheck(m_shape.size() == 3, "");
        DCheck(a < m_shape[0], "");
        DCheck(b < m_shape[1], "");
        DCheck(c < m_shape[2], "");
        return (a * m_shape[1] + b) * m_shape[2] + c;
    }
    size_t offset(size_t a, size_t b, size_t c, size_t d) const {
        DCheck(m_shape.size() == 4, "");
        DCheck(a < m_shape[0], "");
        DCheck(b < m_shape[1], "");
        DCheck(c < m_shape[2], "");
        DCheck(d < m_shape[3], "");
        return ((a * m_shape[1] + b) * m_shape[2] + c) * m_shape[3] + d;
    }

    const T& operator()(tensor_shape index) const { return operator[](offset(index)); }

    // sub-tensor for the fixed value of the first dimension
    Tensor slice(size_t a) {
        DCheck(m_shape.size() > 0, "");
        DCheck(a < m_shape[0], "");
        auto v = m_shape.pop_front().volume();
        return Tensor(m_data + a * v, m_shape.pop_front());
    }

    const Tensor slice(size_t a) const {
        DCheck(m_shape.size() > 0, "");
        DCheck(a < m_shape[0], "");
        auto v = m_shape.pop_front().volume();
        return Tensor(m_data + a * v, m_shape.pop_front());
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

    size_t size() const { return m_vector.size(); }

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
