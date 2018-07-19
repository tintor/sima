#pragma once

static_assert(sizeof(int) == 4);
static_assert(sizeof(long) == 8);
static_assert(sizeof(long long) == 8);

using int128 = __int128;
using uint128 = __uint128_t;

namespace std {

inline int128 abs(int128 a) { return (a < 0) ? -a : a; }

}

template<typename T>
void hash_combine(size_t& seed, const T& value) {
	// The code is from `hash_combine` function of the Boost library. See
	// http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine .
	seed ^= std::hash<T>()(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {

template <typename A, typename B> struct hash<pair<A, B>> {
	size_t operator()(const pair<A, B>& x) const {
		size_t seed = 0;
		hash_combine(seed, x.first);
		hash_combine(seed, x.second);
		return seed;
	}
};

}
