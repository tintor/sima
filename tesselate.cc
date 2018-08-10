#include "tesselate.h"
#include "range.h"
#include "auto.h"
#include "scalar.h"
#include "aabb.h"
#include "is_valid.h"

// D - disjoint
// O - overlap
// V - vertex / vertex touch (could be colinear, but not overlapping)
// X - interior intersection
// A - T intersection: A or B is touching interior of PQ
// B - T intersection: P or Q is touching interior of AB
char relate(ivec2 a, ivec2 b, ivec2 p, ivec2 q) {
	int sp = sign(area(p, a, b));
	int sq = sign(area(q, a, b));
	if (sp == 0 && sq == 0) { // colinear
		if (aabb<ivec2>(a, b).overlaps(aabb<ivec2>(p, q)))
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
bool relate_abxo(ivec2 a, ivec2 b, ivec2 p, ivec2 q, long ab) {
	long pa = edge_area(p, a);
	long bp = edge_area(b, p);
	long sp = pa + ab + bp; // not overflow safe!

	long qa = edge_area(q, a);
	long bq = edge_area(b, q);
	long sq = qa + ab + bq; // not overflow safe!

	if (sp == 0) {
		if (sq == 0) // colinear
			return aabb<ivec2>(a, b).overlaps(aabb<ivec2>(p, q));
		if (equal(a, p) || equal(b, p))
			return false;
	} else {
		if (sq == 0)
			if (equal(a, q) || equal(b, q))
				return false;

		if (sign(sp) * sign(sq) > 0)
			return false;
	}

    long pq = edge_area(p, q);
	long z3 = pq + qa; // not overflow safe!
	int sign_sa = (pa < z3) - (z3 < pa);

	long z4 = bp + pq; // not overflow safe!
	int sign_sb = (bq < z4) - (z4 < bq);
	return sign_sa * sign_sb <= 0;
}

const int EMPTY = std::numeric_limits<int>::min();

bool intersects_polyline(isegment2 s, long ab, ipolygon2::iterator begin, ipolygon2::iterator last) {
	ivec2 p = *last++;
	for (auto it = begin; it != last; it++) {
		ivec2 q = *it;
		if (q.x != EMPTY) {
			assert(!equal(p, q));
			if (relate_abxo(s.a, s.b, p, q, ab))
				return true;
			p = q;
		}
	}
	return false;
}

void tesselate(ipolygon2 poly, imesh2& tess) {
	assert(is_valid(poly));
	auto orig_tess_size = tess.size();
	auto orig_poly_size = poly.size();
	tess.reserve(tess.size() + poly.size() - 2);

	long w = area(poly);
	auto n = poly.size();

	// TODO use bits to avoid modifying poly!

	auto begin = poly.begin();
	auto last = poly.end() - 1;

	auto a = last;
	auto b = begin;
	auto c = begin + 1;
	long ab = edge_area(*a, *b);
	long bc = edge_area(*b, *c);
	long ca = edge_area(*c, *a);
	if (n > 3)
	while (true) {
		long v = ab + bc + ca;
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
