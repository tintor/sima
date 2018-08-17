#include "tesselate.h"
#include "range.h"
#include "auto.h"
#include "scalar.h"
#include "aabb.h"
#include "is_valid.h"

// TODO review entire file after conversion from int -> double

// D - disjoint
// O - overlap
// V - vertex / vertex touch (could be colinear, but not overlapping)
// X - interior intersection
// A - T intersection: A or B is touching interior of PQ
// B - T intersection: P or Q is touching interior of AB
char relate(double2 a, double2 b, double2 p, double2 q) {
	int sp = sign(area(p, a, b));
	int sq = sign(area(q, a, b));
	if (sp == 0 && sq == 0) { // colinear
		if (aabb2(a, b).intersects(aabb2(p, q)))  // TODO overlaps changed to intersects
			return 'O';
		if (equal(a, p) || equal(a, q) || equal(b, p) || equal(b, q))
			return 'V';
		return 'D';
	}

	if (equal(a, p) || equal(a, q) || equal(b, p) || equal(b, q))
		return 'V';

	int sab = sign(area(a, p, q)) * sign(area(b, p, q));
	int spq = sp * sq;
	if (sab < 0 && spq < 0)
		return 'X';
	if (sab == 0 && spq < 0)
		return 'A';
	if (sab < 0 && spq == 0)
		return 'B';
	return 'D';
}

// same as relate(), but faster
bool relate_abxo(double2 a, double2 b, double2 p, double2 q, long ab) {
	double pa = edge_area(p, a);
	double bp = edge_area(b, p);
	double sp = pa + ab + bp;

	double qa = edge_area(q, a);
	double bq = edge_area(b, q);
	double sq = qa + ab + bq;

	if (sp == 0) {
		if (sq == 0) // colinear
			return aabb2(a, b).intersects(aabb2(p, q)); // TODO overlaps changed to intersects
		if (equal(a, p) || equal(b, p))
			return false;
	} else {
		if (sq == 0)
			if (equal(a, q) || equal(b, q))
				return false;

		if (sign(sp) * sign(sq) > 0)
			return false;
	}

    double pq = edge_area(p, q);
    double z3 = pq + qa;
	double z4 = bp + pq;
	int sign_sa = (pa < z3) - (z3 < pa);
	int sign_sb = (bq < z4) - (z4 < bq);
	return sign_sa * sign_sb <= 0;
}

const int EMPTY = std::numeric_limits<int>::min();

bool intersects_polyline(segment2 s, long ab, polygon2::iterator begin, polygon2::iterator last) {
	double2 p = *last++;
	for (auto it = begin; it != last; it++) {
		double2 q = *it;
		if (q.x != EMPTY) {
			assert(!equal(p, q));
			if (relate_abxo(s.a, s.b, p, q, ab))
				return true;
			p = q;
		}
	}
	return false;
}

void tesselate(polygon2 poly, mesh2& tess) {
	assert(is_valid(poly));
	auto orig_tess_size = tess.size();
	auto orig_poly_size = poly.size();
	tess.reserve(tess.size() + poly.size() - 2);

	double w = area(poly);
	auto n = poly.size();

	// TODO use bits to avoid modifying poly!

	auto begin = poly.begin();
	auto last = poly.end() - 1;

	auto a = last;
	auto b = begin;
	auto c = begin + 1;
	double ab = edge_area(*a, *b);
	double bc = edge_area(*b, *c);
	double ca = edge_area(*c, *a);
	if (n > 3)
	while (true) {
		double v = ab + bc + ca;
		if (abs(w) == abs(w - v) + abs(v) && !intersects_polyline({*c, *a}, ca, begin, last)) {
			w -= v;
			tess.emplace_back(*a, *b, *c);
			b->x = EMPTY;
			n -= 1;
			if (b == begin)
				begin = c;
			if (b == last)
				last = a;
			b = c;
			ab = -ca;
			if (n == 3) {
				// find next c
				if (c == last)
					c = begin;
				else {
					c++;
					while (c->x == EMPTY)
						c++;
				}
				break;
			}
			if (n <= (last - begin + 1) - 30) {
				// condense
				auto r = poly.begin(), w = poly.begin();
				while (true) {
					if (r->x != EMPTY) {
						if (w < r)
							*w = *r;
						if (++w == poly.begin() + n)
							break;
					}
					r++;
				}
				begin = poly.begin();
				last = begin + n - 1;
				a = begin;
				b = begin + 1;
				c = begin + 2;
				ab = edge_area(*a, *b);
				bc = edge_area(*b, *c);
				ca = edge_area(*c, *a);
				continue;
			}
		} else {
			a = b;
			b = c;
			ab = bc;
		}
		// find next c
		if (c == last)
			c = begin;
		else {
			c++;
			while (c->x == EMPTY)
				c++;
		}
		bc = edge_area(*b, *c);
		ca = edge_area(*c, *a);
	}
	tess.emplace_back(*a, *b, *c);
	assert(tess.size() == orig_tess_size + orig_poly_size - 2);
}
