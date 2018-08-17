#pragma once
#include "std.h"
#include <functional>
#include <immintrin.h>

template<typename T>
ulong hash(T a) {
	constexpr ulong seed = 0;
	return _mm_crc32_u64(~seed, std::hash<T>()(a));
}

template<typename T, typename U>
ulong hash(T a, U b) {
	return _mm_crc32_u64(std::hash<T>()(a), std::hash<U>()(b));
}

template<typename U>
ulong hash(double a, U b) {
	return _mm_crc32_u64(*reinterpret_cast<ulong*>(&a), std::hash<U>()(b));
}

template<typename T, typename... Args>
ulong hash(T a, Args... args) {
	return _mm_crc32_u64(hash(args...), std::hash<T>()(a));
}

namespace std {

template <typename A, typename B> struct hash<pair<A, B>> {
	size_t operator()(const pair<A, B>& x) const { return hash(x.first, x.second); }
};

}
