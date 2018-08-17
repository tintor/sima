#pragma once
#include "format.h"
#include "variant.h"

// T can be any number-like type
template <typename T> struct range {
	explicit range(T& end) : range(0, end) {}
	explicit range(const T& end) : range(0, end) {}

	range(T begin, T& end, T inc = 1)
		: _begin(begin), _end(0), _end_ptr(&end), _inc(inc) {
		assert((inc > 0 && begin <= end) || (inc < 0 && begin >= end));
	}

	range(T begin, const T& end, T inc = 1)
		: _begin(begin), _end(end), _end_ptr(&_end), _inc(inc) {
		assert((inc > 0 && begin <= end) || (inc < 0 && begin >= end));
	}

	struct iterator {
		iterator(T pos, T inc) : _var(Var0{pos, inc}) {}
		iterator(const T* end_ptr) : _var(end_ptr) {}

		T& operator*() { assert(!is_end()); return pos(); }
		T operator*() const { return is_end() ? end() : pos(); }
		void operator++() { pos() += inc(); }
		void operator++(int) { pos() += inc(); }
		void operator--() { pos() -= inc(); }
		void operator--(int) { pos() -= inc(); }

		bool operator!=(iterator e) {
			if (is_end()) {
				if (e.is_end())
					return end() == e.end();
				return (e.inc() > 0) ? (e.pos() < pos()) : (e.pos() > pos());
			}
			if (e.is_end())
				return (inc() > 0) ? (pos() < e.end()) : (pos() > e.end());
			assert(inc() == e.inc());
			return (inc() > 0) ? (pos() < e.pos()) : (pos() > e.pos());
		}

	  private:
		T& pos() { return mpark::get<0>(_var).pos; }
		T pos() const { return mpark::get<0>(_var).pos; }
		T inc() const { return mpark::get<0>(_var).inc; }
		T end() const { return *mpark::get<1>(_var); }
		bool is_end() const { return _var.index() != 0; }

		struct Var0 {
			T pos, inc;
		};
		mpark::variant<Var0, const T*> _var;
	};

	iterator begin() const { return iterator(_begin, _inc); }
	iterator end() const { return iterator(_end_ptr); }
	T inc() const { return _inc; }

  private:
	const T _begin;
	const T _end;
	const T* _end_ptr;
	const T _inc;
};

template <typename T>
void format_e(string& s, string_view spec, range<T> v) {
	using namespace std::literals;
	s += "range("sv;
	if (*v.begin() != 0 || std::abs(v.inc()) != 1) {
		format_e(s, spec, *v.begin());
		s += ", "sv;
	}
	const auto ei = v.end();
	format_e(s, spec, *ei);
	if (std::abs(v.inc()) != 1) {
		s += ", "sv;
		format_e(s, spec, v.inc());
	}
	s += ')';
}
