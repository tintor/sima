#pragma once

static_assert(sizeof(char) == 1);
static_assert(sizeof(int) == 4);
static_assert(sizeof(long) == 8);
static_assert(sizeof(long long) == 8);

using uchar = unsigned char;

using uint = unsigned int;
using ulong = unsigned long;

using int128 = __int128;
using uint128 = __uint128_t;

using size_t = ulong;

namespace std {

inline int128 abs(int128 a) { return (a < 0) ? -a : a; }

}
