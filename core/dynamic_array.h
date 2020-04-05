#pragma once
#include <core/int.h>
#include <type_traits>
#include <cstdlib>
#include <cstring>
#include <new>
#include <cassert>
#include <initializer_list>

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
