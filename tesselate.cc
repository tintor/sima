#include "tesselate.h"
#include "range.h"
#include "auto.h"
#include "scalar.h"
#include "aabb.h"
#include "is_valid.h"
#include "segment.h"
#include "classify.h"
#include "exception.h"
#include "timestamp.h"

static bool Overlaps(aabb2 a, aabb2 b) {
	double2 mn = vmax(a.min, b.min), mx = vmin(a.max, b.max);
	double2 t{Tolerance, Tolerance};
	return all(mn - mx <= t) && squared(mn - mx) > squared(Tolerance);
}

/*
inline double signed_double_edge_area(double2 a, double2 b) {
	return (a.x + b.x) * (a.y - b.y);
}

inline double signed_double_area(double2 a, double2 b, double2 c) {
	return signed_double_edge_area(a, b) + signed_double_edge_area(b, c) + signed_double_edge_area(c, a);
}

inline int Sign(segment2 s, double2 v) {
	return Sign(signed_double_area(s.a, s.b, v) / length(s.a - s.b));
}
*/

bool fast_relate_abxo(segment2 p, segment2 q) {
	int sp = Sign(p, q.a);
	int sq = Sign(p, q.b);
	if (sp * sq > 0)
		return false;
	if (sp == 0 && sq == 0) // colinear
		return Overlaps(aabb2(p), aabb2(q));
	int sab = Sign(q, p.a) * Sign(q, p.b);
	if (sab > 0)
		return false;
	return sp * sq < 0 || sab < 0;
}

static bool relate_abxo(const polygon2& polygon, segment2 e) {
	for (segment2 s : Edges(polygon))
		if (fast_relate_abxo(e, s))
			return true;
	return false;
}

// 0.253s

template<typename T, typename RND>
T random_int(T lo, T hi, RND& engine) {
	std::uniform_int_distribution<T> dist(lo, hi);
	return dist(engine);
}

vector<uint> make_deck(uint size) {
	std::default_random_engine rnd;
	rnd.seed(0);
	vector<uint> deck(size);
	for (uint i : range(deck.size())) {
		uint j = random_int<uint>(0, i, rnd);
		deck[i] = deck[j];
		deck[j] = i;
	}
	return deck;
}

void compress_deck(vector<uint>& deck, uint cards) {
	uint r = 0;
	uint w = 0;
	while (r < deck.size()) {
		if (deck[r] < cards)
			deck[w++] = deck[r];
		r += 1;
	}
	deck.resize(w);
}

long t_deck = 0;
long t_pip = 0;
long t_relate = 0;
long t_erase = 0;

// slow and simple
void tesselateSimple(const polygon2& in, mesh2& mesh) {
	auto deck = make_deck(in.size());
	uint g = 0;
	uint fails = 0;

	polygon2 polygon = in;
	double a1 = signed_double_area(polygon);
	while (polygon.size() > 3) {
		if (fails == polygon.size()) {
			print("polygon = %s\n", polygon);
			THROW(runtime_error, "no solution found");
		}

		Timestamp ta;
		while (deck[g] >= polygon.size())
			if (++g == deck.size()) {
				compress_deck(deck, polygon.size());
				g = 0;
			}
		uint b = deck[g];
		if (++g == deck.size()) {
			compress_deck(deck, polygon.size());
			g = 0;
		}
		t_deck += ta.elapsed();

		uint a = (b > 0) ? b - 1 : (polygon.size() - 1);
		uint c = (b < polygon.size() - 1) ? b + 1 : 0;
		segment2 e(polygon[a], polygon[c]);

		// O(N)
		Timestamp tb;
		double a2 = signed_double_area(triangle2{e.a, polygon[b], e.b});
		bool pip = (abs(a2) + abs(a1 - a2) - abs(a1)) < 1e-9;
		//bool pip = PointInPolygon((e.a + e.b) / 2, polygon);
		t_pip += tb.elapsed();
		if (!pip) {
			fails += 1;
			continue;
		}
		// O(N)
		Timestamp tc;
		bool rel = relate_abxo(polygon, e);
		t_relate += tc.elapsed();
		if (rel) {
			fails += 1;
			continue;
		}

		fails = 0;
		mesh.emplace_back(e.a, polygon[b], e.b);
		// O(N)
		a1 -= a2;
		Timestamp td;
		polygon.erase(polygon.begin() + b);
		t_erase += td.elapsed();
	}

	assert(polygon.size() == 3);
	mesh.emplace_back(polygon[0], polygon[1], polygon[2]);
}

// TODO review entire file after conversion from int -> double

// same as relate(), but faster
bool relate_abxo(double2 a, double2 b, double2 p, double2 q, double ab) {
	double pa = signed_double_edge_area(p, a);
	double bp = signed_double_edge_area(b, p);
	double sp = pa + ab + bp;

	double qa = signed_double_edge_area(q, a);
	double bq = signed_double_edge_area(b, q);
	double sq = qa + ab + bq;

	// TODO with doubles this will never happen!
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

bool intersects_polyline(segment2 s, double ab, polygon2::iterator begin, polygon2::iterator last) {
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
