#pragma once
#include <core/format.h>
#include <core/property.h>
#include <core/std.h>
#include <core/util.h>

using dim_t = uint;

struct dim4 {
    dim4(dim_t a = 0, dim_t b = 0, dim_t c = 0, dim_t d = 0)
            : dim4{a, b, c, d, ' ', ' ', ' ', ' '} { }

    dim4(dim_t a, dim_t b, dim_t c, dim_t d, char an, char bn, char cn, char dn)
            : d{a, b, c, d}, n{an, bn, cn, dn} {
        Check(a != 0 || (b == 0 && an == ' '));
        Check(b != 0 || (c == 0 && bn == ' '));
        Check(c != 0 || (d == 0 && cn == ' '));
        Check(d != 0 || dn == ' ');

        size = 0;
        elem = 1;
        for (int i = 0; i < 4; i++) if (this->d[i] > 0) {
            elem *= this->d[i];
            size += 1;
        }
    }

    static dim4 empty() {
        dim4 d;
        d.elem = 0;
        return d;
    }

    int ndims() const { return size; }

    dim_t elements() const { return elem; }

    dim_t operator[](uint i) const { Check(i < size); return d[i]; }
    dim_t& operator[](uint i) { Check(i < size); return d[i]; }
    bool operator==(dim4 o) const { return d == o.d; }
    bool operator!=(dim4 o) const { return d != o.d; }

    dim_t back() const { Check(size > 0); return d[size - 1]; }

    dim4 pop_front() const { return {d[1], d[2], d[3], 0, n[1], n[2], n[3], ' '}; }

    dim4 pop_back() const {
        int r = rank();
        if (r <= 0) return *this;
        dim4 e = *this;
        e[r - 1] = 0;
        return e;
    }

    dim4 push_front(dim_t a) const { Check(d[3] == 0); return {a, d[0], d[1], d[2], ' ', n[0], n[1], n[2]}; }

    // 0 scalar : all dimensions are 1
    // 1 vector : if one dimension is >1 and rest are 1
    // 2 matrix
    // 3 cuboid
    // 4 hyper-cuboid
    int rank() const {
        int r = 0;
        for (int i = 0; i < 4; i++) if (d[i] > 1) r += 1;
        return r;
    }

    string str() const {
        string s;
        s += '[';
        for (int i = 0; i < size; i++) {
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
    dim_t elem;
    uint size;
};

// Tensor is a pointer to multi-dimensional array (doesn't own memory)
template <typename T = float>
class Tensor {
   public:
    using type = T;

    Tensor() : m_data(nullptr), m_shape(dim4::empty()) {}
    Tensor(T* data, dim4 shape) : m_data(data), m_shape(shape) {
        Check(data || shape.elements() == 0);
        Check(shape.elements() != 0 || data == nullptr);
    }

    Tensor(const Tensor<T>& o) : m_data(o.m_data), m_shape(o.m_shape) {}

    operator bool() const { return m_data != nullptr; }

    T* data() { return m_data; }
    T const* data() const { return m_data; }

    T* begin() { return data(); }
    T* end() { return data() + size; }
    const T* begin() const { return data(); }
    const T* end() const { return data() + size; }

    Property(Tensor) {
        operator size_t() const { return PropertyT<Tensor, __LINE__ - 1>::parent->m_shape.elements(); }
    }
    size;

    const auto& shape() const { return m_shape; }

    Property(Tensor) {
        operator uint() const { return PropertyT<Tensor, __LINE__ - 1>::parent->m_shape.rank(); }
    }
    rank;

    T& operator[](size_t index) {
        DCheck(index < size(), format("%s < %s", index, size));
        return m_data[index];
    }
    const T& operator[](size_t index) const {
        DCheck(index < size(), format("%s < %s", index, size));
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
        return m_shape == o.m_shape && std::equal(m_data, m_data + size, o.m_data);
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
        for (size_t i = 0; i < size; i++) {
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
        m_vector.resize(o.size);
        Tensor<T>::m_data = m_vector.data();
        // TODO not safe if strides are added to Tensor
        std::copy(o.data(), o.data() + o.size, m_vector.data());
        return *this;
    }

    TProperty(Size, VTensor) {
        // using Prop = PropertyT<VTensor, __LINE__ - 1>;
        operator size_t() const { return Size::parent->m_vector.size(); }
    }
    size;

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
