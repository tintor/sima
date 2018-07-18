#include "tesselate.h"
#include "range.h"
#include "catch.hpp"

inline long m(int a, int b, int c, int d) {
	return (long(a) + b) * (long(c) - d);
}

static int sign(const ipolygon2& poly) {
	long area = 0;
    auto a = poly.back();
	for (auto b : poly) {
        area += m(a.x, b.x, a.y, b.y);
		a = b;
	}
	return (0 < area) - (area < 0);
}

static int sign(ivec2 a, ivec2 b, ivec2 c) {
	long area = 0;
    area += m(a.x, b.x, a.y, b.y);
    area += m(b.x, c.x, b.y, c.y);
    area += m(c.x, a.x, c.y, a.y);
	return (0 < area) - (area < 0);
}

imesh2 tesselate(ipolygon2 poly) {
	REQUIRE(poly.size() >= 3);
	imesh2 m;
	m.reserve(poly.size() - 2);
	auto sign_poly = sign(poly);
	REQUIRE(sign_poly != 0);
	ipolygon2 copy;
	copy.reserve(poly.size() - 1);
	while (poly.size() > 2) {
		auto a = poly[poly.size() - 2];
		auto b = poly.back();
		for (auto c : poly) {
			if (c == poly.back())
				c = copy[0];
			if (sign_poly == sign(a, b, c)) {
				m.emplace_back(a, b, c);
			} else {
				copy.push_back(b);
				if (poly.size() - copy.size() < 3)
					return m;
				a = b;
			}
			b = c;
		}
		REQUIRE(copy.size() < poly.size());
		poly.clear();
		std::swap(poly, copy);
	}
	return m;
}
