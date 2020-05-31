#pragma once
#include <core/int.h>
#include <core/std.h>
#include <immintrin.h>

#include <functional>

// Hash combine using a single SSE4.2 instruction!
struct hash {
    ulong seed = ~static_cast<ulong>(0);
    hash& operator<<(ulong b) {
        seed = _mm_crc32_u64(seed, b);
        return *this;
    }
};

inline hash operator<<(hash h, int a) {
    ulong b = a;
    return h << b;
}

inline hash operator<<(hash h, float a) {
    ulong b = *reinterpret_cast<uint*>(&a);
    return h << b;
}

inline hash operator<<(hash h, double a) { return h << *reinterpret_cast<ulong*>(&a); }

template <typename A, typename B>
hash operator<<(hash h, const pair<A, B>& v) {
    return h << v.first << v.second;
}

template <typename T>
struct hash_t {
    size_t operator()(const T& a) const { return (hash() << a).seed; }
};

/*namespace std {

template<typename T>
struct hash {
        size_t operator()(const T& x) const {
                return (::hash() << x).seed;
        }
};

}*/
