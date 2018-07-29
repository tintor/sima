#pragma once
#include <vector>

template<typename A, typename B>
class zip {
public:
	zip(std::vector<A>& a, std::vector<B>& b) {
		_a = a.data();
		_b = b.data();
		_n = std::min(a.size(), b.size());
	}

	struct iter {
		A* _a;
		B* _b;
		bool operator!=(iter e) { return _a != e._a || _b != e._b; }
		void operator++() { _a++; _b++; }
		std::pair<A&, B&> operator*() { return {*_a, *_b}; }
	};

	iter begin() { return {_a, _b}; }
	iter end() { return {_a + _n, _b + _n}; }

private:
	A* _a;
	B* _b;
	size_t _n;
};

template<typename A, typename B>
class czip {
public:
	czip(const std::vector<A>& a, const std::vector<B>& b) {
		_a = a.data();
		_b = b.data();
		_n = std::min(a.size(), b.size());
	}

	struct iter {
		const A* _a;
		const B* _b;
		bool operator!=(iter e) { return _a != e._a || _b != e._b; }
		void operator++() { _a++; _b++; }
		std::pair<A, B> operator*() { return {*_a, *_b}; }
	};

	iter begin() const { return {_a, _b}; }
	iter end() const { return {_a + _n, _b + _n}; }

private:
	const A* _a;
	const B* _b;
	size_t _n;
};

