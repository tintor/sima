#pragma once

#include <core/std.h>
#include <core/range.h>

template<typename Fn, typename FnPrime>
inline double NewtonMethod(double x, int iters, const Fn& f, const FnPrime& f_prime) {
	for (auto i : range(iters))
		x -= f(x) / f_prime(x);
	return x;
}

// TODO Newton method for x = vector

