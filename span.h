#pragma once
#include "int.h"
#include <iterator>
#include <vector>
#include <array>
#include <cassert>

// View into an array or container of contiguous elements,
// with optional stride allowing skipping elements and reversing order.
// Allows function to operate on "array" of elements without caring about memory ownership.
//
// sizeof(span) == 16 bytes
template<typename T>
class span {
public:
	using M = typename std::remove_const<T>::type;

	constexpr span(T* data, size_t size, ssize_t stride = 1)
			: _data(data), _size(size), _stride(stride) {
		assert(stride != 0);
		assert(size <= std::numeric_limits<uint>::max());
		assert(stride <= std::numeric_limits<int>::max());
		assert(stride >= std::numeric_limits<int>::min());
	}
	constexpr span(T* begin, T* end) : span(begin, end - begin) { }
	constexpr span(const span& a) : _data(a._data), _size(a._size), _stride(a._stride) { }

	template<size_t N>
	constexpr span(const std::array<M, N>& v) : span(v.data(), v.data() + v.size()) { }
	constexpr span(const std::vector<M>& v) : span(v.data(), v.data() + v.size()) { }

	constexpr T& operator[](size_t idx) const { assert(idx < size()); return _data[idx * stride()]; }

	constexpr T* data() const { return _data; }
	constexpr size_t size() const { return _size; }
	constexpr ssize_t stride() const { return _stride; }

	constexpr bool empty() const { return size() == 0; }

	constexpr T* begin() const { return _data; }
	constexpr T* end() const { return _data + size() * stride(); }

	constexpr T& front() { return begin()[0]; }
	constexpr T& back() { return end()[-1]; }

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

	constexpr span first(size_t count) const { return subspan(0, count); }
	constexpr span last(size_t count) const { assert(count <= size()); return subspan(size() - count, count); }
	constexpr span reverse() const { return subspan(size() - 1, size(), -stride()); }

private:
	T* _data;
	uint _size;
	int _stride;
};

/*template<typename T, size_t N>
span<const T> Span(const std::array<T, N>& v) { return {v.data(), v.data() + N}; }

template<typename T>
span<const T> Span(const std::vector<T>& v) { return {v.data(), v.data() + v.size()}; }
*/
