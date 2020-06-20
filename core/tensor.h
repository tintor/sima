#pragma once
#include <core/format.h>
#include <core/property.h>
#include <core/std.h>
#include <core/util.h>

using dim_t = uint;

struct dim4 {
    dim4(dim_t a = 0, dim_t b = 0, dim_t c = 0, dim_t d = 0, char an = ' ', char bn = ' ', char cn = ' ', char dn = ' ')
            : d{a, b, c, d}, n{an, bn, cn, dn} {
        Check(a != 0 || (b == 0 && an == ' '));
        Check(b != 0 || (c == 0 && bn == ' '));
        Check(c != 0 || (d == 0 && cn == ' '));
        Check(d != 0 || dn == ' ');
    }

    int ndims() const {
        if (d[3] != 0) return 4;
        if (d[2] != 0) return 3;
        if (d[1] != 0) return 2;
        if (d[0] != 0) return 1;
        return 0;
    }

    size_t elements() const {
        if (d[3] != 0) return size_t(d[0]) * size_t(d[1]) * size_t(d[2]) * size_t(d[3]);
        if (d[2] != 0) return size_t(d[0]) * size_t(d[1]) * size_t(d[2]);
        if (d[1] != 0) return size_t(d[0]) * size_t(d[1]);
        if (d[0] != 0) return d[0];
        return 1;
    }

    char name(uint i) const { Check(i < 4 && n[i] != 0); return n[i]; }
    dim_t operator[](uint i) const { Check(i < 4 && n[i] != 0); return d[i]; }
    bool operator==(dim4 o) const { return d == o.d && n == o.n; }
    bool operator!=(dim4 o) const { return d != o.d || n != o.n; }

    dim_t back() const { Check(d[0] != 0); return d[ndims() - 1]; }

    dim4 pop_front() const { return {d[1], d[2], d[3], 0, n[1], n[2], n[3], ' '}; }

    dim4 pop_back() const {
        if (d[0] == 0) return *this;
        dim4 e = *this;
        auto size = ndims();
        e.d[size - 1] = 0;
        e.n[size - 1] = ' ';
        return e;
    }

    dim4 push_front(dim_t a, char an = ' ') const {
        Check(d[3] == 0);
        return {a, d[0], d[1], d[2], an, n[0], n[1], n[2]};
    }

    dim4 push_back(dim_t a, char an = ' ') const {
        Check(d[3] == 0);
        Check(a > 0);
        dim4 e = *this;
        auto size = ndims();
        e.d[size] = a;
        e.n[size] = an;
        return e;
    }

    string str() const {
        string s;
        s += '[';
        for (int i = 0; i < 4 && d[i] != 0; i++) {
            if (i != 0) s += ' ';
            if (n[i] != ' ') s += n[i];
            format_s(s, "%s", d[i]);
        }
        s += ']';
        return s;
    }

    operator string() const { return str(); }

private:
    std::array<dim_t, 4> d;
    std::array<char, 4> n;
};

// Tensor is a pointer to multi-dimensional array (doesn't own memory)
template <typename T = float>
class Tensor {
   public:
    using type = T;

    Tensor() : m_data(nullptr), m_shape() {}
    Tensor(T* data, dim4 shape) : m_data(data), m_shape(shape) {
        Check(data || shape.ndims() == 0);
    }

    Tensor(const Tensor<T>& o) : m_data(o.m_data), m_shape(o.m_shape) {}

    operator bool() const { return m_data != nullptr; }

    T* data() { return m_data; }
    T const* data() const { return m_data; }

    T* begin() { return data(); }
    T* end() { return data() + size(); }
    const T* begin() const { return data(); }
    const T* end() const { return data() + size(); }

    auto ndims() const { return m_shape.ndims(); }
    auto size() const { return m_data ? m_shape.elements() : 0; }

    dim_t dim(uint i) const { return m_shape[i]; }
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

    const T& operator()(dim4 index) const { return operator[](offset(index)); }

    // sub-tensor for the fixed value of the first dimension
    Tensor slice(size_t a) {
        DCheck(m_shape.size() > 0, "");
        DCheck(a < m_shape[0], "");
        size_t v = m_shape.pop_front().elements();
        return Tensor(m_data + a * v, m_shape.pop_front());
    }

    const Tensor slice(size_t a) const {
        DCheck(m_shape.size() > 0, "");
        DCheck(a < m_shape[0], "");
        size_t v = m_shape.pop_front().elements();
        return Tensor(m_data + a * v, m_shape.pop_front());
    }

    void copy_from(const Tensor& o) {
        if (shape() != o.shape()) Fail(format("%s vs %s", string(shape()), string(o.shape())));
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
    dim4 m_shape;
};

// VTensor is multi-dimensional vector (owns memory)
template <typename T>
class VTensor : public Tensor<T> {
   public:
    VTensor() {}
    VTensor(dim4 shape, T init = T()) : m_vector(shape.elements(), init) {
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
        o.m_shape = dim4();
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

    auto size() const { return m_vector.size(); }

    void reshape(dim4 new_shape, T init = T()) {
        Tensor<T>::m_shape = new_shape;
        m_vector.resize(new_shape.elements(), init);
        Tensor<T>::m_data = m_vector.data();
    }

   private:
    vector<T> m_vector;
};

using tensor = Tensor<float>;
using vtensor = VTensor<float>;
