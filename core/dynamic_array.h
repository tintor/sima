#pragma once
#include <core/int.h>
#include <type_traits>
#include <cstdlib>
#include <new>
#include <cassert>

template<typename T, typename = std::enable_if_t<std::is_scalar<T>::value || std::is_pod<T>::value>>
class raw_array {
public:
	raw_array() : _data(nullptr) { }

	raw_array(uint size, bool zero_init = false) {
		_data = reinterpret_cast<T*>(zero_init ? calloc(sizeof(T), size) : malloc(sizeof(T) * size));
		if (_data == nullptr)
			throw std::bad_alloc();
	}

	~raw_array() { free(_data); }

	raw_array(raw_array&& e) {
		_data = e._data;
		e._data = nullptr;
	}

	raw_array& operator=(raw_array&& e) {
		if (this != &e) {
			free(_data);
			_data = e._data;
			e._data = nullptr;
		}
		return *this;
	}

	raw_array(const raw_array& e) = delete;

	raw_array& operator=(const raw_array& e) = delete;

	T& operator[](uint idx) {
		return _data[idx];
	}

	const T& operator[](uint idx) const {
		return _data[idx];
	}

	void resize(uint size) {
		if (size == 0) {
			free(_data);
			_data = nullptr;
			return;
		}
		T* data = reinterpret_cast<T*>(realloc(_data, sizeof(T) * size));
		if (!data)
			throw std::bad_alloc();
		_data = data;
	}

	T* data() const { return _data; }

private:
    T *_data;
};

template<typename T, typename = std::enable_if_t<std::is_scalar<T>::value> >
class dynamic_array {
public:
	dynamic_array() : _data(nullptr), _size(0) { }

	dynamic_array(uint size, bool zero_init = false) {
		_data = reinterpret_cast<T*>(zero_init ? calloc(sizeof(T), size) : malloc(sizeof(T) * size));
		if (_data == nullptr)
			throw std::bad_alloc();
		_size = size;
	}

	dynamic_array(dynamic_array&& e) {
		_data = e._data;
		_size = e._size;
		e._data = nullptr;
		e._size = 0;
	}

	dynamic_array(const dynamic_array& e) : _data(nullptr), _size(0) {
		resize(e.size());
		memcpy(_data, e.begin(), _size * sizeof(T));
	}

	~dynamic_array() { free(_data); }

	dynamic_array& operator=(const dynamic_array& e) {
		if (this != &e) {
			resize(e.size());
			memcpy(_data, e.begin(), _size * sizeof(T));
		}
		return *this;
	}

	dynamic_array& operator=(dynamic_array&& e) {
		if (this != &e) {
			_data = e._data;
			_size = e._size;
			e._data = nullptr;
			e._size = 0;
		}
		return *this;
	}

	T& operator[](uint idx) {
		assert(idx < _size);
		return _data[idx];
	}

	T operator[](uint idx) const {
		assert(idx < _size);
		return _data[idx];
	}

	void resize(uint size) {
		if (size == _size)
			return;
		if (size == 0) {
			free(_data);
			_data = nullptr;
			_size = 0;
			return;
		}
		T* data = reinterpret_cast<T*>(realloc(_data, sizeof(T) * size));
		if (!data)
			throw std::bad_alloc();
		_data = data;
		_size = size;
	}

	void swap(dynamic_array& o) {
		std::swap(_data, o._data);
		std::swap(_size, o._size);
	}

	uint size() const { return _size; }
	T* begin() { return _data; }
	T* end() { return _data + _size; }
private:
    T *_data;
    uint _size;
};
