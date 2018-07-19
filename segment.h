#pragma once
#include "glm.h"

template<typename V>
struct segment {
	V a, b;

	V& operator[](int idx) { return (&a)[idx]; }
	V operator[](int idx) const { return (&a)[idx]; }

	bool operator==(const segment& v) const { return a == v.a && b == v.b; }
	bool operator!=(const segment& v) const { return a != v.a || b != v.b; }

	segment reversed() const { return {b, a}; }
};

template<typename T>
void format_e(std::string& s, std::string_view spec, segment<T> v) {
	format_e(s, "", v.a);
	s += ", ";
	format_e(s, "", v.b);
}

using isegment2 = segment<ivec2>;
using isegment3 = segment<ivec3>;

using segment2_ref = segment<ivec2&>;
using segment3_ref = segment<ivec3&>;

namespace std {

template<typename V> struct hash<segment<V>> {
	size_t operator()(const segment<V>& s) const {
		size_t seed = 0;
		hash_combine(seed, s.a);
		hash_combine(seed, s.b);
		return seed;
	}
};

}
