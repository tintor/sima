#pragma once
#include <vector>
#include <array>

// TODO add stride (to be able to skip elements or go in reverse)
// TODO add sub_array method
template<typename T>
struct array_ptr {
	T* _begin;
	T* _end;

	array_ptr(T* b, T* e) {
		_begin = b;
		_end = e;
		auto bi = reinterpret_cast<size_t>(b);
		auto ei = reinterpret_cast<size_t>(e);
		assert((bi == 0 && ei == 0) || (bi != 0 && ei != 0 && (ei - bi) % sizeof(T) == 0));
	}

	template<size_t N>
	array_ptr(std::array<T, N>& v) {
		_begin = v.data();
		_end = v.data() + v.size();
	}

	array_ptr(std::vector<T>& v) {
		_begin = v.data();
		_end = v.data() + v.size();
	}

	T& operator[](size_t idx) {
		return _begin[idx];
	}
	
	const T& operator[](size_t idx) const {
		return _begin[idx];
	}

	size_t size() const {
		return _end - _begin;
	}
	
	T* begin() {
		return _begin;
	}
	
	T* end() {
		return _end;
	}
};

template<typename T>
struct array_cptr {
	const T* _begin;
	const T* _end;

	array_cptr(const T* b, const T* e) {
		_begin = b;
		_end = e;
	}
	
	array_cptr(array_ptr<T> a) {
		_begin = a._begin;
		_end = a._end;
	}

	template<size_t N>
	array_cptr(const std::array<T, N>& v) {
		_begin = v.data();
		_end = v.data() + v.size();
	}

	array_cptr(const std::vector<T>& v) {
		_begin = v.data();
		_end = v.data() + v.size();
	}

	const T& operator[](size_t idx) const {
		return _begin[idx];
	}

	size_t size() const {
		return _end - _begin;
	}

	const T* begin() const {
		return _begin;
	}
	
	const T* end() const {
		return _end;
	}
};

