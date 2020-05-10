#pragma once
#include <core/int.h>
#include <core/vector.h>

#include <array>
#include <cassert>
#include <iterator>

// View into an array or container of contiguous elements,
// with optional stride allowing skipping elements and reversing order.
// Allows function to operate on "array" of elements without caring about memory ownership.
//
// sizeof(span) == 16 bytes
template <typename T>
class span {
   public:
    using M = typename std::remove_const<T>::type;
    using iterator = T*;
    using const_iterator = const T*;

    operator span<const T>() { return span<const T>(data(), size()); }

    constexpr span(T* data, size_t size, ssize_t stride = 1) : _data(data), _size(size), _stride(stride) {
        assert(stride != 0);
        assert(size <= std::numeric_limits<uint>::max());
        assert(stride <= std::numeric_limits<int>::max());
        assert(stride >= std::numeric_limits<int>::min());
    }
    constexpr span(T* begin, T* end) : span(begin, end - begin) {}
    constexpr span(const span& a) : _data(a._data), _size(a._size), _stride(a._stride) {}

    template <size_t N>
    constexpr span(const std::array<M, N>& v) : span(v.data(), v.data() + v.size()) {}

    template <size_t N>
    constexpr span(std::array<T, N>& v) : span(v.data(), v.data() + v.size()) {}

    template <typename A>
    constexpr span(const vector<M, A>& v) : span(v.data(), v.data() + v.size()) {}

    template <typename A>
    constexpr span(vector<T, A>& v) : span(v.data(), v.data() + v.size()) {}

    constexpr span(const std::initializer_list<M>& v) : span(v.begin(), v.end()) {}

    constexpr T& operator[](size_t idx) const {
        assert(idx < size());
        return _data[idx * stride()];
    }

    constexpr bool operator==(span v) const noexcept {
        if (size() != v.size()) return false;
        for (size_t i = 0; i < size(); i++)
            if (operator[](i) != v[i]) return false;
        return true;
    }
    constexpr bool operator!=(span v) const noexcept { return !operator==(v); }

    constexpr T* data() const { return _data; }
    constexpr size_t size() const { return _size; }
    constexpr ssize_t stride() const { return _stride; }

    constexpr bool empty() const { return size() == 0; }

    constexpr T* begin() const { return _data; }
    constexpr T* end() const { return _data + size() * stride(); }

    constexpr T& front() { return begin()[0]; }
    constexpr T& back() { return end()[-1]; }

    constexpr const T& front() const { return begin()[0]; }
    constexpr const T& back() const { return end()[-1]; }

    constexpr span subspan(size_t begin, size_t new_size, ssize_t new_stride = 1) const {
        assert(new_stride != 0);
        if (new_size > 0) {
            if (new_stride >= 0)
                assert(begin + new_size * new_stride <= size());
            else
                assert(begin >= (new_size - 1) * -new_stride);
        }
        return span(_data + begin * stride(), new_size, stride() * new_stride);
    }

    constexpr span pop_front(size_t n = 1) const { return subspan(n, size() - n); }
    constexpr span pop_back(size_t n = 1) const { return subspan(0, size() - n); }

    constexpr span first(size_t count) const { return subspan(0, count); }
    constexpr span last(size_t count) const {
        assert(count <= size());
        return subspan(size() - count, count);
    }
    constexpr span reverse() const { return subspan(size() - 1, size(), -stride()); }

   private:
    T* _data;
    uint _size;
    int _stride;
};

template <typename T>
using cspan = span<const T>;
