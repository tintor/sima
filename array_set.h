#pragma once
#include "dynamic_array.h"
#include "bits.h"
#include "range.h"
#include "catch.hpp"
#include <functional>
#include <iostream>

// T must have these methods:
// no destructor
// bool empty()
// void clear()
// operator==(K key)
// std::hash

template<typename T>
class open_addressing_hashtable {
public:
	open_addressing_hashtable() { }
	open_addressing_hashtable(const open_addressing_hashtable&) = delete;
	open_addressing_hashtable& operator=(const open_addressing_hashtable&) = delete;

	template<typename Key>
	T* get(Key a) {
		assert(consistent());
		uint i = index(a);
		while (!empty(i)) {
			if (equals(i, a))
				return &_array[i];
			i = (i + 1) & _mask;
		}
		return nullptr;
	}
	
	template<typename Key>
	const T* get(Key a) const {
		assert(consistent());
		uint i = index(a);
		while (!empty(i)) {
			if (equals(i, a))
				return &_array[i];
			i = (i + 1) & _mask;
		}
		return nullptr;
	}

	// insert or update
	// returns true on insert, false on update
	bool insert_or_update(T a) {
		assert(consistent());
		if (_size >= capacity()) {
			rehash(std::max<uint>(8, (_mask + 1) * 2) - 1);
			assert(consistent());
		}
		uint i = index(a);
		while (!empty(i)) {
			if (equals(i, a)) {
				assert(consistent());
				return false;
			}
			i = (i + 1) & _mask;
		}
		_array[i] = a;
		_size += 1;
		assert(consistent());
		return true;
	}

	// never inserts
	bool update(T a) {
		assert(consistent());
		uint i = index(a);
		while (!empty(i)) {
			if (equals(i, a)) {
				_array[i] = a;
				assert(consistent());
				return true;
			}
			i = (i + 1) & _mask;
		}
		assert(consistent());
		return false;
	}

	// never updates
	bool insert(T a) {
		throw "not implemented";
	}

	void insert_unsafe(T a) {
		assert(!contains(a));
		assert(_array.data() != nullptr);
		assert(_size < capacity());
		uint i = index(a);
		while (empty(i))
			i = (i + 1) & _mask;
		_array[i] = a;
		_size += 1;
		assert(consistent());
	}

	template<typename Key>
	bool erase(Key a) {
		assert(consistent());
		uint i = index(a);
		uint z = i;
		while (!empty(i)) {
			if (equals(i, a)) {
				_size -= 1;
				fill_gap(z, i);
				assert(consistent());
				return true;
			}
			i = (i + 1) & _mask;
		}
		assert(consistent());
		return false;
	}
	
	template<typename Key>
	void erase_unsafe(Key a) {
		assert(get(a) != nullptr);
		uint i = index(a);
		uint z = i;
		while (true) {
			if (equals(i, a)) {
				_size -= 1;
				fill_gap(z, i);
				assert(consistent());
				return;
			}
			i = (i + 1) & _mask;
		}
	}

	void reserve(uint size) {
		assert(consistent());
		uint mask = round_up_power2(size * 2) - 1;
		if (mask + 1 > _mask + 1)
			rehash(mask);
		assert(consistent());
	}

	void shrink() {
		assert(consistent());
		if (_size == 0) {
			_array.resize(0);
			_mask = -1;
			assert(consistent());
			return;
		}
		uint mask = std::max(8, round_up_power2(_size * 2)) - 1;
		if (mask + 1 < _mask + 1)
			rehash(mask);	
		assert(consistent());
	}

	void clear() {
		assert(consistent());
		if (_size > 0) {
			for (auto i : range(_mask + 1))
				clear(i);
			_size = 0;
		}
		assert(consistent());
	}

	uint size() const { return _size; }
	uint capacity() const { return (_mask + 1) / 2; }

	class iterator {
	public:
		iterator(T* pos) : _pos(pos) { }
		T& operator*() { assert(!_pos->empty()); return *_pos; }
		void operator++() {
			_pos += 1;
			while (_pos->empty())
				_pos += 1;
		}
		bool operator!=(iterator e) const { return _pos != e._pos; }
	private:
		T* _pos;
	};

	iterator begin() {
		assert(consistent());
		if (_size == 0)
			return iterator(nullptr);
		T* d = _array.data();
		while (d->empty())
			d += 1;
		return iterator(d);
	}

	iterator end() {
		assert(consistent());
		if (_size == 0)
			return iterator(nullptr);
		T* e = _array.data() + _mask;
		while (e->empty())
			e -= 1;
		return iterator(e + 1);
	}

	// DEBUG ONLY
	void print() const {
		std::cout << "[";
		for (auto i : range(_mask + 1))
			std::cout << " " << _array[i];
		std::cout << "]" << std::endl;
	}
	
	uint max_probes() const {
		/*uint a = 0;
		if (_array[0] != empty) {
			while (_array[(a - 1) & _mask] != empty)
				a -= 1;
		}
		for (uint i = a; i != (a - 1) & _mask; i = (i + 1) & _mask) {
			if (_array
		}*/
		return 0;
	}

private:
	bool empty(uint i) const { return _array[i].empty(); }
	
	template<typename Key>
	bool equals(uint i, Key a) const { return _array[i].equals(a); }
	
	void clear(uint i) { _array[i].clear(); }

	static ulong fmix64(ulong k) {
		k ^= k >> 33;
		k *= 0xff51afd7ed558ccdllu;
		k ^= k >> 33;
		k *= 0xc4ceb9fe1a85ec53llu;
		k ^= k >> 33;
		return k;
	}

	// TODO take high order bits instead of low order bits so that on regrow we get space between every two slots
	template<typename Key>
	static uint index(Key a, uint mask) {
		ulong h = fmix64(std::hash<Key>()(a));
		uint p = h;
		uint q = h >> 32;
		return (p ^ q) & mask;
	}
	
	template<typename Key>
	uint index(Key a) const {
		return index(a, _mask);
	}

	bool internal_contains(T a) const {
		assert(!a.empty());
		uint i = index(a);
		while (!empty(i)) {
			if (equals(i, a))
				return true;
			i = (i + 1) & _mask;
		}
		return false;
	}

	bool consistent() const {
		// trivial checks
		REQUIRE((bool)(_mask == -1 || _array.data() != nullptr));
		REQUIRE((bool)(_mask != -1 || _array.data() == nullptr));
		REQUIRE(_size <= (_mask + 1) / 2);
		// _size must match number of non-empty slots
		uint size = 0;
		for (auto i : range(_mask + 1))
			if (!empty(i))
				size += 1;
		REQUIRE(size == _size);
		// every element must be in proper slot
		for (auto i : range(_mask + 1))
			if (!empty(i))
				assert(internal_contains(_array[i]));
		// no duplicates allowed
		if (_size > 1) {
			raw_array<T> copy(_size);
			uint w = 0;
			for (auto i : range(_mask + 1))
				if (!empty(i))
					copy[w++] = _array[i];
			std::sort(copy.data(), copy.data() + _size);
			for (auto i : range(_size - 1))
				assert(!copy[i].equals(copy[i + 1]));
		}
		return true;
	}

	static bool between(uint a, uint b, uint c) {
		if (a <= b)
			return b <= c || c <= a;
		return b <= c && c <= a;
	}

	void fill_gap(uint a, uint b) {
		uint c = (b + 1) & _mask;
		if (empty(c)) {
			clear(b);
			return;
		}

		// decrease a to find start of block
		while (!empty((a - 1) & _mask))
			a = (a - 1) & _mask;

		// find element right from b that can be moved to slot b
		while (true) {
			uint e = index(_array[c]);
			if (between(a, e, b)) {
				_array[b] = _array[c];
				b = c;
			}
			c = (c + 1) & _mask;
			if (empty(c))
				break;
		}
		clear(b);
	}

	void rehash(uint mask) {
		assert(is_power2(mask + 1));
		// allocate and init new array
		raw_array<T> array(mask + 1);
		for (auto i : range(mask + 1))
			clear(i);
		// copy from _array to array
		for (auto j : range(_mask + 1)) {
			uint i = index(_array[j], mask);
			while (!empty(i))
				i = (i + 1) & mask;
			array[i] = _array[j];
		}
		// commit new array
		_array = std::move(array);
		_mask = mask;
	}

private:
	raw_array<T> _array;
	uint _size = 0;
	uint _mask = -1;
};

template<typename T, typename Value, T empty_key = 0>
class array_map {
private:	
	struct Slot {
		T key;
		Value value;
		bool equals(T k) const { return key == k; }
		bool equals(Slot s) const { return key == s.key; }
		bool empty() { return key == empty_key; }
		void clear() { key = empty_key; }
	};

public:
	array_map() { }
	array_map(const array_map&) = delete;
	array_map& operator=(const array_map&) = delete;

	bool contains(T a) const {
		return _map.contains(a);
	}

	bool insert(T a, Value v) {
		return _map.insert({a, v});
	}
	
	void insert_unsafe(T a, Value v) {
		return _map.insert_unsafe({a, v});
	}

	bool erase(T a) {
		return _map.erase(a);
	}
	
	void erase_unsafe(T a) {
		return _map.erase_unsafe(a);
	}

	void reserve(uint size) {
		_map.reserve(size);
	}

	void shrink() {
		_map.shrink();
	}

	void clear() {
		_map.clear();
	}

	uint size() const { return _map.size(); }
	uint capacity() const { return _map.capacity(); }
	
	class iterator {
	public:
		iterator(T* data) : _data(data) { }
		T operator*() { assert(!_data->empty()); return *_data; }
		void operator++() {
			_data += 1;
			while (_data->empty())
				_data += 1;
		}
		bool operator!=(iterator e) const { return _data != e._data; }
	private:
		array_map::Slot* _data;
	};

	iterator begin() {
		return iterator(_map._begin());
	}

	iterator end() {
		return iterator(_map._end());
	}

	uint max_probes() {
		return _map.max_probes();
	}
	
	void print() {
		_map.print();
	}

private:
	open_addressing_hashtable<Slot> _map;
};

namespace std {

template<typename T, typename Value>
struct hash<typename ::array_map<T, Value>::Slot> {
	size_t operator()(const array_map<T, Value>::Slot& s) const {
		return hash<T>()(s.key);
	}
};

}

template<typename T, typename Value>
std::ostream& operator<<(std::ostream&& os, const typename array_map<T, Value>::Slot& s) {
	return os << s.key << ":" << s.value;
}

template<typename T, T empty_value = 0>
class array_set {
private:
	struct Slot {
		T key;
		bool equals(T k) const { return key == k; }
		bool equals(Slot s) const { return key == s.key; }
		bool empty() { return key == empty_value; }
		void clear() { key = empty_value; }
	};

public:
	array_set() { }
	array_set(const array_set&) = delete;
	array_set& operator=(const array_set&) = delete;

	bool contains(T a) const {
		return _map.get(a) != nullptr;
	}

	bool insert(T a) {
		return _map.insert(Slot{a});
	}
	
	void insert_unsafe(T a) {
		return _map.insert_unsafe(Slot{a});
	}

	bool erase(T a) {
		return _map.erase(a);
	}
	
	void erase_unsafe(T a) {
		return _map.erase_unsafe(a);
	}

	void reserve(uint size) {
		_map.reserve(size);
	}

	void shrink() {
		_map.shrink();
	}

	void clear() {
		_map.clear();
	}

	uint size() const { return _map.size(); }
	uint capacity() const { return _map.capacity(); }

	class iterator {
	public:
		iterator(T* data) : _data(data) { }
		T operator*() { assert(_data->empty()); return *_data; }
		void operator++() {
			_data += 1;
			while (_data->empty())
				_data += 1;
		}
		bool operator!=(iterator e) const { return _data != e._data; }
	private:
		Slot* _data;
	};

	iterator begin() {
		return iterator(_map._begin());
	}

	iterator end() {
		return iterator(_map._end());
	}

	uint max_probes() {
		return _map.max_probes();
	}

	void print() {
		_map.print();
	}

private:
	open_addressing_hashtable<Slot> _map;
};

template<typename T, T empty>
std::ostream& operator<<(std::ostream&& os, typename array_set<T, empty>::Slot s) {
	return os << s.key;
}
