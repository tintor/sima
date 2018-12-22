#pragma once
struct nullopt_t {
};

constexpr nullopt_t nullopt = nullopt_t();

template<typename T>
class optional {
public:
	optional() : _has_value(false) { }
	optional(nullopt_t) : _has_value(false) { }
	optional(T value) : _has_value(true), _value(std::move(value)) { }
	~optional() {
		if (_has_value)
			_value.~T();
	}
	optional(const optional& o) : _has_value(o._has_value) {
		if (_has_value) {
			new(&_value) T(o._value);
		}
	}
	optional(optional&& o) : _has_value(o._has_value) {
		if (_has_value) {
			new(&_value) T(std::move(o._value));
		}
		o._has_value = false;
	}
	void operator=(const optional& o) {
		if (_has_value && o._has_value) {
			_value.operator=(o._value);
			return;
		}
		if (_has_value) {
			_value.~T();
			_has_value = false;
			return;
		}
		if (o._has_value) {
			new(&_value) T(o._value);
			_has_value = true;
			return;
		}
	}
	void operator=(optional&& o) {
		if (_has_value && o._has_value) {
			_value.operator=(std::move(o._value));
			o._has_value = false;
			return;
		}
		if (_has_value) {
			_value.~T();
			_has_value = false;
			return;
		}
		if (o._has_value) {
			new(&_value) T(std::move(o._value));
			_has_value = true;
			o._has_value = false;
			return;
		}
	}
	void operator=(nullopt_t) {
		if (_has_value) {
			_value.~T();
			_has_value = false;
		}
	}
	void operator=(T value) {
		if (_has_value) {
			_value.operator=(std::move(value));
			return;
		}
		new(&_value) T(std::move(value));
		_has_value = true;
	}
	bool has_value() const { return _has_value; }
	T& operator*() { assert(_has_value); return _value; }
	T* operator->() { assert(_has_value); return &_value; }
private:
	union {
		T _value;
		char _dummy;
	};
	bool _has_value;
};
