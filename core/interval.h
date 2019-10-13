#pragma once
#include <core/util.h>
#include <core/std.h>

template<typename T = double>
struct interval {
	T min, max;

	interval() {
		min = std::numeric_limits<double>::max();
		max = std::numeric_limits<double>::min();
	}

	interval(T min, T max) : min(min), max(max) {
	}

	void add(T a) {
		minimize(min, a);
		maximize(max, a);
	}

	interval operator+(double a) const {
		return interval(min + a, max + a);
	}

	interval operator-(double a) const {
		return interval(min - a, max - a);
	}
};
