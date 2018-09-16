#include "catch.hpp"
#include "classify.h"

TEST_CASE("Sign(segment2, double2)", "[classify]") {
	segment2 s(double2{2, 0}, double2{1, 0});
	REQUIRE(Sign(s, double2{1.5, Tolerance}) == 1);
	REQUIRE(Sign(s, double2{1.5, Tolerance * 0.999999}) == 0);

	REQUIRE(Sign(s, double2{1.5, 0.1}) == 1);
	REQUIRE(Sign(s, double2{1.5, 0}) == 0);
	REQUIRE(Sign(s, double2{1.5, -0.1}) == -1);
	REQUIRE(Sign(s, double2{1, 0}) == 0);
	REQUIRE(Sign(s, double2{2, 0}) == 0);
	REQUIRE(Sign(s, double2{0, 0}) == 0);

	segment2 p(double2{3, 2}, double2{3, 1});
	REQUIRE(Sign(p, double2{3 + Tolerance, 1.5}) == -1);
	REQUIRE(Sign(p, double2{3 + Tolerance * 0.999999, 1.5}) == 0);

	REQUIRE(Sign(p, double2{3.1, 1.5}) == -1);
	REQUIRE(Sign(p, double2{3, 1.5}) == 0);
	REQUIRE(Sign(p, double2{-2.9, 1.5}) == 1);
	REQUIRE(Sign(p, double2{3, 1}) == 0);
	REQUIRE(Sign(p, double2{3, 2}) == 0);
	REQUIRE(Sign(p, double2{3, 0}) == 0);

	std::default_random_engine rnd;
	for (int i = 0; i < 10000; i++) {
		double2 a = uniform2(rnd, -10, 10);
		double2 b = uniform2(rnd, -10, 10);
		double2 c = uniform2(rnd, -10, 10);
		segment2 s(a, b);
		REQUIRE(Sign(s, a) == 0);
		REQUIRE(Sign(s, b) == 0);
		REQUIRE(Sign(s.reversed(), c) == -Sign(s, c));
		for (int j = 0; j < 100; j++)
			REQUIRE(Sign(s, s.linear(uniform(rnd, -2, 3))) == 0);
	}
}

TEST_CASE("operator<<(vector<T>, span<T>)", "[classify]") {
	vector<int> a;
	a << span<const int>{2, 3};
	a << span<const int>{};
	a << span<const int>{4};
	REQUIRE(a == vector<int>{2, 3, 4});
}

double2 swap(double2 a) {
	return double2{a.y, a.x};
}

// TODO test ray/edge crossing logic
// TODO test "shape" of Tolerance area around the border, is it round or square? (sample many points and estimate area of Tolerance area)
TEST_CASE("Classify(xpolygon2, double2)", "[classify]") {
	xpolygon2 a;
	a.add({double2{0, 0}, double2{1, 0}, double2{1, 1}, double2{0, 1}});

	REQUIRE(Classify(a, double2{1, 1}) == 0);
	REQUIRE(Classify(a, double2{1.1, 1.1}) == 1);
	REQUIRE(Classify(a, double2{0.9, 1}) == 0);
	REQUIRE(Classify(a, double2{0.9, 0.9}) == -1);

	std::default_random_engine rnd;
	for (int i = 0; i < 10000; i++) {
		print("i=%s\n", i);
		double2 a = uniform2(rnd, -10, 10);
		double2 b = uniform2(rnd, -10, 10);
		double2 c = uniform2(rnd, -10, 10);
		xpolygon2 poly;
		poly.add(a);
		poly.add(b);
		poly.add(c);

		xpolygon2 poly2;
		poly2.add(c);
		poly2.add(b);
		poly2.add(a);

		xpolygon2 poly3;
		poly3.add(-a);
		poly3.add(-b);
		poly3.add(-c);

		xpolygon2 poly4;
		poly4.add(swap(a));
		poly4.add(swap(b));
		poly4.add(swap(c));

		for (double2 p : {a, b, c}) {
			print("[%s %s %s] %s\n", a, b, c, p);
			REQUIRE(Classify(poly, p) == 0);
		   	REQUIRE(Classify(poly2, p) == 0);
			REQUIRE(Classify(poly3, -p) == 0);
			REQUIRE(Classify(poly4, swap(p)) == 0);
		}
		for (int j = 0; j < 100; j++) {
			double t = uniform(rnd, -1, 2);
			for (double2 p : {
					segment2(a, b).linear(t),
					segment2(b, c).linear(t),
					segment2(c, a).linear(t),
					uniform2(rnd, -10, 10)}) {
				int k = Classify(poly, p);
				REQUIRE(k == Classify(poly2, p));
				REQUIRE(k == Classify(poly3, -p));
				REQUIRE(k == Classify(poly4, swap(p)));
			}

			double s = uniform(rnd, 0, 1);
			double2 p = segment2(a, b).linear(s);
			REQUIRE(Classify(poly, p) == 0);
			REQUIRE(Classify(poly2, p) == 0);
			REQUIRE(Classify(poly3, -p) == 0);
			REQUIRE(Classify(poly4, swap(p)) == 0);

			s = uniform(rnd, 0.001, 0.999);
			p = segment2(p, c).linear(s);
			REQUIRE(Classify(poly, p) == -1);
			REQUIRE(Classify(poly2, p) == -1);
			REQUIRE(Classify(poly3, -p) == -1);
			REQUIRE(Classify(poly4, swap(p)) == -1);
		}
	}
}

template<typename T>
bool L_L(T a, T b, T c) {
	return a < b && b < c;
}

template<typename T>
bool LE_LE(T a, T b, T c) {
	return a <= b && b <= c;
}

template<typename T>
bool LE_L(T a, T b, T c) {
	return a <= b && b < c;
}

template<typename T>
bool L_LE(T a, T b, T c) {
	return a < b && b <= c;
}

template<typename T>
bool E_E(T a, T b, T c) {
	return a == b && b == c;
}

TEST_CASE("relate(segment2, segment2)", "[classify]") {
	segment2 a(double2{1, 0}, double2{1, 1});
	segment2 b(double2{1, 0.1}, double2{1, 0.9});
	segment2 c(double2{2, 0}, double2{2, 1});
	REQUIRE(relate(a, b) == 'O');
	REQUIRE(relate(a, c) == 'D');

	auto verify = [](double2 a, double2 b, double2 c, double2 d) {
		char k = relate(segment2(a, b), segment2(c, d));
		REQUIRE(k == relate(segment2(b, a), segment2(c, d)));
		REQUIRE(k == relate(segment2(b, a), segment2(d, c)));
		REQUIRE(k == relate(segment2(a, b), segment2(d, c)));

		// TODO swap x,y
		// TODO flip x
		// TODO flip y
		return k;
	};

	std::default_random_engine rnd;
	for (int i = 0; i < 1000000; i++) {
		double2 a = uniform2(rnd, -10, 10);
		double2 b = uniform2(rnd, -10, 10);
		while (Equals(a, b))
			b = uniform2(rnd, -10, 10);
		double2 c = uniform2(rnd, -10, 10);
		while (Equals(a, c) || Equals(b, c))
			c = uniform2(rnd, -10, 10);
		double2 d = uniform2(rnd, -10, 10);
		while (Equals(a, d) || Equals(b, d) || Equals(c, d))
			d = uniform2(rnd, -10, 10);

		double t = uniform(rnd, -1, 2);
		double2 ab = segment2(a, b).linear(t);
		while (Equals(a, ab) || Equals(b, ab) || Equals(c, ab) || Equals(d, ab)) {
			t = uniform(rnd, -1, 2);
			ab = segment2(a, b).linear(t);
		}

		double s = uniform(rnd, -1, 2);
		double2 ab2 = segment2(a, b).linear(s);
		while (Equals(a, ab2) || Equals(b, ab2) || Equals(c, ab2) || Equals(d, ab2) || Equals(ab, ab2)) {
			s = uniform(rnd, -1, 2);
			ab2 = segment2(a, b).linear(s);
		}

		// TODO two parallel (non-colinear) segments

		// 4 random points
		verify(a, b, c, d);

		// 3 random points
		char k = verify(a, b, b, c);
		REQUIRE(bool(k == 'V' || k == 'O'));

		// 2 random points
		REQUIRE(verify(a, b, a, b) == 'O');

		// 3 random points + 1 random colinear
		verify(a, b, c, ab);

		// 2 random points + 1 random colinear
		k = verify(a, b, b, ab);
		if (t > 1 + Tolerance)
			REQUIRE(k == 'V');
		else
			REQUIRE(k == 'O');

		// 2 random colinear segments
		k = verify(a, b, ab, ab2);
		if (s > 1 + Tolerance && t > 1 + Tolerance)
			REQUIRE(k == 'D');
		else if (s < -Tolerance && t < -Tolerance)
			REQUIRE(k == 'D');
		else if (L_L(Tolerance, s, 1 - Tolerance) || L_L(Tolerance, t, 1 - Tolerance))
			REQUIRE(k == 'O');
		else
			REQUIRE(bool(k == 'O' || k == 'V'));
	}
}

TEST_CASE("Classify(xpolygon2, segment2)", "[classify]") {
	xpolygon2 a;
	a.add({double2{0, 0}, double2{1, 0}, double2{1, 1}, double2{0, 1}});

	REQUIRE(Classify(a, segment2(double2{2, 1}, double2{2, 0})) == 1);

	REQUIRE(Classify(a, segment2(double2{1, 1}, double2{2, 1})) == 0);
	REQUIRE(Classify(a, segment2(double2{1, 1}, double2{1, 0})) == 0);
	REQUIRE(Classify(a, segment2(double2{1, 0.1}, double2{1, 0.9})) == 0);
	REQUIRE(Classify(a, segment2(double2{1, 0.1}, double2{1, 1})) == 0);
	REQUIRE(Classify(a, segment2(double2{0.5, 1}, double2{1.5, 1})) == 0);
	REQUIRE(Classify(a, segment2(double2{-1, 1}, double2{2, 1})) == 0);

	REQUIRE(Classify(a, segment2(double2{1, 1}, double2{0, 0})) == -1);
	REQUIRE(Classify(a, segment2(double2{1, 1}, double2{0.5, 0.5})) == -1);
	REQUIRE(Classify(a, segment2(double2{1.5, 1.5}, double2{0.5, 0.5})) == -1);
	REQUIRE(Classify(a, segment2(double2{0.5, 0.5}, double2{1.5, 0.5})) == -1);
	REQUIRE(Classify(a, segment2(double2{0.5, 0.5}, double2{1, 0.5})) == -1);
	REQUIRE(Classify(a, segment2(double2{0.1, 0.5}, double2{0.9, 0.5})) == -1);
}

TEST_CASE("Classify(xpolygon2, xpolygon2)", "[classify]") {
	xpolygon2 a;
	a.add({double2{0, 0}, double2{1, 0}, double2{1, 1}, double2{0, 1}});
	for (int i = -1; i <= 1; i++) {
		xpolygon2 b;
		double e = 1 + i * 0.1;
		b.add({double2{e, 0}, double2{e+1, 0}, double2{e+1, 1}, double2{e, 1}});
		REQUIRE(Classify(a, b) == i);
	}
}
