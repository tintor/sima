#include "catch.hpp"
#include "classify.h"
#include "primitives.h"

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

double distanceRingPoint(span<const double2> ring, double2 p) {
	double dist = std::numeric_limits<double>::max();
	for (const segment2& e : Edges(ring))
		dist = min(dist, distance(e, p));
	return dist;
}

TEST_CASE("Classify(xpolygon2, double2) simple", "[classify]") {
	xpolygon2 a;
	a.add({double2{0, 0}, double2{1, 0}, double2{1, 1}, double2{0, 1}});

	REQUIRE(Classify(a, double2{1, 1}) == 0);
	REQUIRE(Classify(a, double2{1.1, 1.1}) == 1);
	REQUIRE(Classify(a, double2{0.9, 1}) == 0);
	REQUIRE(Classify(a, double2{0.9, 0.9}) == -1);
}

TEST_CASE("Classify(xpolygon2, double2) random triangle", "[classify]") {
	std::default_random_engine rnd;
	for (int i = 0; i < 10000; i++) {
		double2 a = uniform2(rnd, -10, 10);
		double2 b = uniform2(rnd, -10, 10);
		while (Equals(a, b))
			b = uniform2(rnd, -10, 10);
		double2 c = uniform2(rnd, -10, 10);
		while (Equals(a, c) || Equals(b, c))
			c = uniform2(rnd, -10, 10);

		xpolygon2 poly;
		poly.add({a, b, c});

		xpolygon2 poly2;
		poly2.add({c, b, a});

		xpolygon2 poly3;
		poly3.add({-a, -b, -c});

		xpolygon2 poly4;
		poly4.add({swap(a), swap(b), swap(c)});

		vector<pair<double2, int>> tests;
		tests.emplace_back(a, 0);
		tests.emplace_back(b, 0);
		tests.emplace_back(c, 0);

		for (int j = 0; j < 100; j++) {
			tests.emplace_back(uniform2(rnd, -10, 10), 9);

			double e = uniform(rnd, -1, -0.01);
			tests.emplace_back(segment2(a, b).linear(e), 1);
			tests.emplace_back(segment2(b, a).linear(e), 1);

			double t = uniform(rnd, 0.01, 0.99);
			tests.emplace_back(segment2(a, b).linear(t), 0);

			double s = uniform(rnd, 0.01, 0.99);
			tests.emplace_back(segment2(segment2(a, b).linear(t), c).linear(s), -1);
		}

		for (const auto& [p, cc] : tests) {
			const int c = (cc == 9) ? Classify(poly, p) : cc;
			REQUIRE(Classify(poly, p) == c);
		   	REQUIRE(Classify(poly2, p) == c);
			REQUIRE(Classify(poly3, -p) == c);
			REQUIRE(Classify(poly4, swap(p)) == c);
		}
	}
}

long squared(int2 a) {
	return (long)a.x * a.x + (long) a.y * a.y;
}

template<typename T>
void tsp_sort(vector<T>& poly) {
	// order points into polygon using traveling salesman heuristic
	bool done = false;
	auto size = poly.size();
	while (!done) {
		done = true;
		for (size_t j : range<size_t>(1, size))
			for (size_t i : range(j - 1)) {
				auto c = poly[j];
				auto d = poly[(j + 1) % size];
				auto a = poly[i];
				auto b = poly[(i + 1) % size];
				if (squared(a - b) + squared(c - d) > squared(a - c) + squared(b - d)) {
					// TODO faster if we replace polygon with XOR-linked list
					std::reverse(poly.begin() + i + 1, poly.begin() + j + 1);
					done = false;
				}
			}
	}
}

template<typename RNG>
static vector<int2> random_int_polygon(int size, RNG& rng) {
	unordered_set<int2, hash_t<int2>, equal_t<int2>> points;
	while (points.size() < size) {
		std::uniform_int_distribution uniform(0, 9);
		points.insert(int2{uniform(rng), uniform(rng)});
	}

	vector<int2> poly(points.begin(), points.end());
	tsp_sort(poly);
	return poly;
}


template<typename RNG>
static polygon2 random_polygon(int size, RNG& rng) {
	polygon2 poly;
	while (poly.size() < size) {
		auto p = uniform2(rng, -10, 10);
		if (All(poly, [p](double2 a){return !Equals(a, p);}))
			poly.push_back(p);
	}
	tsp_sort(poly);
	return poly;
}

double2 transform(double2 p, int t) {
	if (t & 1) p.x = -p.x;
	if (t & 2) p.y = -p.y;
	if (t & 4) p = {p.y, p.x};
	return p;
}

void transform(vector<double2>& poly, int t) {
	for (double2& p : poly)
		p = transform(p, t);
	if (t & 8)
		for (int i = 0; i < poly.size() / 2; i++)
			swap(poly[i], poly[poly.size() - 1 - i]);
}

TEST_CASE("Classify(polygon2, double2) random small int polygon", "[classify]") {
	std::default_random_engine rnd;
	polygon2 poly[16];
	for (int i = 0; i < 10000; i++) {
		poly[0].clear();
		for (int2 p : random_int_polygon(8, rnd))
			poly[0].push_back(double2{static_cast<double>(p.x), static_cast<double>(p.y)});
		for (int i = 1; i < 16; i++) {
			poly[i] = poly[0];
			transform(poly[i], i);
		}
		for (int x = 0; x <= 9; x++)
			for (int y = 0; y <= 9; y++) {
				double2 p{static_cast<double>(x), static_cast<double>(y)};
				int cn[3] = {0, 0, 0};
				for (int t = 0; t < 16; t++) {
					double2 tp = transform(p, t);
					cn[1 + Classify(poly[t], tp)] += 1;
				}
				int m = max(cn[0], cn[1], cn[2]);
				int c;
				if (cn[0] == m) c = -1;
				if (cn[1] == m) c = 0;
				if (cn[2] == m) c = 1;
				if (m != 16) {
					print("outside %s, border %s, inside %s\n", cn[0], cn[1], cn[2]);
					for (int t = 0; t < 16; t++) {
						double2 tp = p;
						transform(tp, t);
						print("[%s] Classify(%s, %s) = %s\n", t, "x"/*poly[t]*/, tp, Classify(poly[t], tp));
					}
					REQUIRE(m == 16);
				}
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
		for (int t = 0; t < 8; t++) {
			double2 A = a, B = b, C = c, D = d;
			transform(A, t);
			transform(B, t);
			transform(C, t);
			transform(D, t);
			REQUIRE(k == relate(segment2(A, B), segment2(C, D)));
			REQUIRE(k == relate(segment2(B, A), segment2(C, D)));
			REQUIRE(k == relate(segment2(B, A), segment2(D, C)));
			REQUIRE(k == relate(segment2(A, B), segment2(D, C)));
		}
		return k;
	};

	std::default_random_engine rnd;
	for (int i = 0; i < 300000; i++) {
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

TEST_CASE("Classify(xpolygon2, segment2) square", "[classify]") {
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

TEST_CASE("Classify(xpolygon2, segment2) concave", "[classify]") {
	xpolygon2 a;
	a.add({double2{0, 0}, double2{1, 1}, double2{2, 0}, double2{2, 2}, double2{0, 2}});

	REQUIRE(Classify(a, segment2(double2{0, 0}, double2{2, 0})) == 0);
}

TEST_CASE("Classify(polygon2, segment2) random polyogn", "[classify]") {
	std::default_random_engine rnd;
	polygon2 poly[16];
	for (int i = 0; i < 500; i++) {
	   	poly[0]	= random_polygon(8, rnd);
		for (int t = 1; t < 16; t++) {
			poly[t] = poly[0];
			transform(poly[t], t);
		}

		const auto& p = poly[0];
		vector<pair<segment2, int>> tests;
		tests.emplace_back(segment2(p[0], p[1]), 0);
		for (auto j : range<size_t>(1, p.size() - 1))
			tests.emplace_back(segment2(p[0], p[j]), 9);

		for (int j = 0; j < 100; j++) {
			double2 a = uniform2(rnd, -10, 10);
			double2 b = uniform2(rnd, -10, 10);
			while (Equals(a, b))
				b = uniform2(rnd, -10, 10);

			tests.emplace_back(segment2(a, b), 9);
			tests.emplace_back(segment2(a, p[0]), 9);
			tests.emplace_back(segment2(a, double2{p[0].y, b.y}), 9);
			tests.emplace_back(segment2(double2{a.x, p[0].y}, double2{b.x, p[0].y}), 9);

			double2 c = segment2(p[0], p[1]).linear(uniform(rnd, -0.5, 1.5));
			tests.emplace_back(segment2(a, c), 9);
		}

		for (const auto& [s, cc] : tests) {
			int c = (cc == 9) ? Classify(poly[0], s) : cc;
			for (int t = 0; t < 16; t++) {
				segment2 ts(transform(s.a, t), transform(s.b, t));
				REQUIRE(c == Classify(poly[t], ts));
			}
		}
	}
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

double2 Rotate(double2 v, double a) {
	double c = cos(a), s = sin(a);
	return { v.x * c + v.y * s, v.y * c - v.x * s };
}

void Rotate(polygon2& m, double a) {
	for (double2& v : m)
		v = Rotate(v, a);
}

void Translate(polygon2& m, double2 t) {
	for (double2& v : m)
		v += t;
}

void Scale(polygon2& m, double2 s) {
	for (double2& v : m)
		v *= s;
}

polygon2 MakeRect(double x, double y) {
	return polygon2{double2{x, y}, double2{-x, y}, double2{-x, -y}, double2{x, -y}};
}

polygon2 MakeTriangle(double x, double y) {
	return polygon2{double2{0, y}, double2{-x/2, 0}, double2{x/2, 0}};
}

polygon2 MakeU(double x, double y, double a) {
	return polygon2{double2{-x, -y}, double2{x, -y}, double2{x, y}, double2{x-a, y},
			        double2{x-a, a-y}, double2{a-x, a-y}, double2{a-x, y}, double2{-x, y}};
}

polygon2 MakeL(double x, double y, double a) {
	return polygon2{double2{-x, -y}, double2{x, -y}, double2{x, a-y}, double2{a-x, a-y},
				    double2{a-x, y}, double2{-x, y}};
}

double TransformUntilContact(polygon2& a, const polygon2& b, std::function<void(polygon2&, double)> transform, double inc) {
	double t = 0;
	polygon2 orig = a;
	while (true) {
		a = orig;
		transform(a, t);
		int k = Classify(a, b);
		if (k == 0)
			return t;
		if (k == -1)
			break;
		t += inc;
	}
	// solution is between t-inc and t
	double lower = t - inc, upper = t;
	while (true) {
		REQUIRE(upper - lower > 1e-8);
		t = (lower + upper) / 2;
		a = orig;
		transform(a, t);
		int k = Classify(a, b);
		if (k == 0)
			return t;
		if (k == -1)
			upper = t;
		else
			lower = t;
	}
}

double RotateUntilContact(polygon2& a, const polygon2& b, double t) {
	return TransformUntilContact(a, b, [t](polygon2& m, double a) { Rotate(m, t); }, M_PI / 72);
}

double TranslateUntilContact(polygon2& a, const polygon2& b, double2 t) {
	return TransformUntilContact(a, b, [t](polygon2& m, double a) { Translate(m, t * a); }, 1);
}

double ScaleUntilContact(polygon2& a, const polygon2& b, double t) {
	return TransformUntilContact(a, b, [t](polygon2& m, double a) { Scale(m, 1 + t); }, 0.001);
}

static bool Overlaps(aabb2 a, aabb2 b) {
	double2 mn = vmax(a.min, b.min), mx = vmin(a.max, b.max);
	return all(mn <= mx) && squared(mx - mn) > squared(Tolerance);
}

TEST_CASE("Overlaps", "[classify2]") {
	segment2 a(double2{3.14, 1}, double2{3.14, -1});
	segment2 b(double2{3.14, -2}, double2{3.14, 2});
	REQUIRE(relate(a, b) == 'O');
	double2 t;
	REQUIRE(relate(a, b, nullptr, &t) == 'O');
	REQUIRE(t.x == 0.25);
	REQUIRE(t.y == 0.75);
}

TEST_CASE("Classify(polygon2, polygon2) small box, large box", "[classify2]") {
	polygon2 a = MakeRect(-1, 1);
	polygon2 b = MakeRect(-2, 2);
	Translate(b, double2{5.14, 0});
	double t = TranslateUntilContact(a, b, double2{1, 0});
	vector<Contact2> contacts;
	Classify(a, b, &contacts);
	for (const Contact2& c : contacts)
		print("sa %s, sb %s, normal %s\n", c.sa, c.sb, c.normal);
}

TEST_CASE("Classify(polygon2, polygon2) offset boxes", "[classify2]") {
	polygon2 a = MakeRect(-1, 1);
	polygon2 b = MakeRect(-1, 1);
	Translate(b, double2{5.14, 1});
	double t = TranslateUntilContact(a, b, double2{1, 0});
	vector<Contact2> contacts;
	Classify(a, b, &contacts);
	for (const Contact2& c : contacts)
		print("sa %s, sb %s, normal %s\n", c.sa, c.sb, c.normal);
}

TEST_CASE("Classify(polygon2, polygon2) rotated box", "[classify2]") {
	polygon2 a = MakeRect(-1, 1);
	Rotate(a, M_PI / 4);
	polygon2 b = MakeRect(-1, 1);
	Translate(b, double2{5.14, 0});
	double t = TranslateUntilContact(a, b, double2{1, 0});
	vector<Contact2> contacts;
	Classify(a, b, &contacts);
	for (const Contact2& c : contacts)
		print("sa %s, sb %s, normal %s\n", c.sa, c.sb, c.normal);
}

TEST_CASE("Classify(polygon2, polygon2) same boxes", "[classify2]") {
	polygon2 a = MakeRect(-1, 1);
	polygon2 b = MakeRect(-1, 1);
	Translate(b, double2{5.14, 0});
	double t = TranslateUntilContact(a, b, double2{1, 0});
	vector<Contact2> contacts;
	Classify(a, b, &contacts);
	for (const Contact2& c : contacts)
		print("sa %s, sb %s, normal %s\n", c.sa, c.sb, c.normal);
}

TEST_CASE("Classify(polygon2, polygon2) square and L", "[classify2]") {
	polygon2 a = MakeRect(-1, 1);
	polygon2 b = MakeL(1, 1, 0.2);
	Rotate(b, M_PI * 3 / 2);
	Translate(b, double2{5, 0});
	double t = TranslateUntilContact(a, b, double2{1, 0});
	vector<Contact2> contacts;
	Classify(a, b, &contacts);
	for (const Contact2& c : contacts)
		print("sa %s, sb %s, normal %s\n", c.sa, c.sb, c.normal);
}
