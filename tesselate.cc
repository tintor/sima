#include "tesselate.h"
#include "range.h"
#include "auto.h"
#include "scalar.h"
#include "aabb.h"
#include "is_valid.h"
#include "segment.h"
#include "classify.h"

static bool relate_abxo(const xpolygon2& polygon, segment2 e) {
	for (segment2 s : Edges(polygon))
		if (relate_abxo(e, s))
			return true;
	return false;
}

// slow and simple
void tesselate(const xpolygon2& in, mesh2& mesh) {
	std::default_random_engine rnd;
	rnd.seed(0);

	xpolygon2 polygon = in;
	while (polygon.vertices().size() > 3) {
		std::uniform_int_distribution<size_t> dist(0, polygon.vertices().size() - 1);
		size_t b = dist(rnd);
		uint r = polygon.ring_of_vertex(b);

		uint rel = b - polygon.offset(r);
		uint size = polygon[r].size();
		uint a = (rel > 0) ? b - 1 : (polygon.offset(r) + size - 1);
		uint c = (rel < size - 1) ? b + 1 : polygon.offset(r);

		segment2 e(polygon.vertices()[a], polygon.vertices()[c]);

		// O(N)
		if (relate_abxo(polygon, e))
			continue;
		// O(N)
		if (Classify(polygon, (e.a + e.b) / 2, /*check_edges*/false) >= 0)
			continue;

		mesh.emplace_back(e.a, polygon.vertices()[b], e.b);
		// O(N)
		polygon.remove_vertex(b);
	}

	const auto& v = polygon.vertices();
	assert(polygon.size() == 1);
	assert(v.size() == 3);
	mesh.emplace_back(v[0], v[1], v[2]);
}

// TODO review entire file after conversion from int -> double

static bool Overlaps(aabb2 a, aabb2 b) {
	double2 mn = vmin(a.max, b.max), mx = vmax(a.min, b.min);
	return all(mn <= mx) && squared(mx - mn) > squared(Tolerance);
}

// same as relate(), but faster
bool relate_abxo(double2 a, double2 b, double2 p, double2 q, long ab) {
	double pa = signed_double_edge_area(p, a);
	double bp = signed_double_edge_area(b, p);
	double sp = pa + ab + bp;

	double qa = signed_double_edge_area(q, a);
	double bq = signed_double_edge_area(b, q);
	double sq = qa + ab + bq;

	if (sp == 0) {
		if (sq == 0) // colinear
			return Overlaps(aabb2(a, b), aabb2(p, q));
		if (equal(a, p) || equal(b, p))
			return false;
	} else {
		if (sq == 0)
			if (equal(a, q) || equal(b, q))
				return false;

		if (sign(sp) * sign(sq) > 0)
			return false;
	}

    double pq = signed_double_edge_area(p, q);
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

	double w = signed_double_area(poly);
	auto n = poly.size();

	// TODO use bits to avoid modifying poly!

	auto begin = poly.begin();
	auto last = poly.end() - 1;

	auto a = last;
	auto b = begin;
	auto c = begin + 1;
	double ab = signed_double_edge_area(*a, *b);
	double bc = signed_double_edge_area(*b, *c);
	double ca = signed_double_edge_area(*c, *a);
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
				ab = signed_double_edge_area(*a, *b);
				bc = signed_double_edge_area(*b, *c);
				ca = signed_double_edge_area(*c, *a);
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
		bc = signed_double_edge_area(*b, *c);
		ca = signed_double_edge_area(*c, *a);
	}
	tess.emplace_back(*a, *b, *c);
	assert(tess.size() == orig_tess_size + orig_poly_size - 2);
}
