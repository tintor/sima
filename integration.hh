template<typename T>
T EulersMethod(T x0, real t0, real dt, std::function<T(T x, real t)> derivative) {
	T k1 = dt * derivative(x0, t0);
	return x0 + k1;
}

template<typename T>
T MidpointMethod(T x0, real t0, real dt, std::function<T(T x, real t)> derivative) {
	T k1 = dt * derivative(x0, t0);
	T k2 = dt * derivative(x0 + k1 / 2, t0 + dt / 2);
	return x0 + k2;
}

// computes X[dt] given X[0] and derivative function of X at x0 and t
// T can be scalar, vector or matrix
template<typename T>
T RungeKutta4(T x0, real t0, real dt, std::function<T(T x, real d)> derivative) {
	T k1 = dt * derivative(x0, t0);
	T k2 = dt * derivative(x0 + k1 / 2, t0 + dt / 2);
	T k3 = dt * derivative(x0 + k2 / 2, t0 + dt / 2);
	T k4 = dt * derivative(x0 + k3, t0 + dt);
	return x0 + (k1 + 2 * (k2 + k3) + k4) / 6;
}

TEST_CASE("Numerical integrators") {
	// differential edquation: y' = y, y(0) = 1
	// exact solution: y(t) = e^t
	// simulating 320 steps from 0.0 to 4.0, step size 0.0125
	const auto& der = [](real e, real t) { return e; };
	real k = 4.0 / 320;

	real e = 1;
	FOR(i, 320) e = EulersMethod<real>(e, i * k, k, der);
	CHECK(abs(e - exp(4.0)) <= 1.34);

	e = 1;
	FOR(i, 320) e = MidpointMethod<real>(e, i * k, k, der);
	CHECK(abs(e - exp(4.0)) <= 0.0057);

	e = 1;
	FOR(i, 320) e = RungeKutta4<real>(e, i * k, k, der);
	CHECK(abs(e - exp(4.0)) <= 4.4e-08);
}
