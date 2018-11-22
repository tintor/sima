#pragma once
#include <functional>

// T = vector, S = scalar
template<typename T, typename S>
T EulersMethod(T x0, S t0, S dt, std::function<T(T x, S t)> derivative) {
	T k1 = dt * derivative(x0, t0);
	return x0 + k1;
}

template<typename T, typename S>
T MidpointMethod(T x0, S t0, S dt, std::function<T(T x, S t)> derivative) {
	T k1 = dt * derivative(x0, t0);
	T k2 = dt * derivative(x0 + k1 / 2., t0 + dt / 2.);
	return x0 + k2;
}

// computes X[dt] given X[0] and derivative function of X at x0 and t
// T can be scalar, vector or matrix
template<typename T, typename S>
T RungeKutta4(T x0, S t0, S dt, std::function<T(T x, S d)> derivative) {
	T k1 = dt * derivative(x0, t0);
	T k2 = dt * derivative(x0 + k1 / 2., t0 + dt / 2.);
	T k3 = dt * derivative(x0 + k2 / 2., t0 + dt / 2.);
	T k4 = dt * derivative(x0 + k3, t0 + dt);
	return x0 + (k1 + 2. * (k2 + k3) + k4) / 6.;
}

// For solving systems of variable size, without memory allocation.
void BigRungeKutta4(const double* x0, double t0, double dt, double* k1, double* k2, double* k3, double* k4, double* out, uint count,
		std::function<void(const double* x, double d, double* d, int count)> derivative) {
	derivative(x0, t0, k1);

	for (uint i = 0; i < count; i++)
		out[i] = x0[i] + k1[i] * (dt / 2);
	derivative(out, t0 + dt / 2, k2);

	for (uint i = 0; i < count; i++)
		out[i] = x0[i] + k2[i] * (dt / 2);
	derivative(out, t0 + dt / 2, k3);

	for (uint i = 0; i < count; i++)
		out[i] = x0[i] + k3[i] * dt;
	derivative(out, t0 + dt, k4);

	for (uint i = 0; i < count; i++)
		out[i] = x0[i] + (k1[i] + k4[i]) * (dt / 6) + (k2[i] + k3[i]) * (dt / 3);
}
