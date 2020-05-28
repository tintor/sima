#pragma once
#include <core/format.h>
#include <core/range.h>
#include <arrayfire.h>

// -1 empty : if any dimension is 0
// 0 scalar : all dimensions are 1
// 1 vector : if one dimension is >1 and rest are 1
// 2 matrix
// 3 cuboid
// 4 hyper-cuboid
inline int rank(af::dim4 a) {
    int r = 0;
    for (int i = 0; i < 4; i++) {
        if (a[i] == 0) return -1;
        if (a[i] > 1) r += 1;
    }
    return r;
}

inline string to_string(af::dim4 m) {
    ulong a = m[0], b = m[1], c = m[2], d = m[3];
    if (d != 1) return format("[%s %s %s %s]", a, b, c, d);
    if (c != 1) return format("[%s %s %s]", a, b, c);
    if (b != 1) return format("[%s %s]", a, b);
    if (a != 1) return format("[%s]", a);
    return "[]";
}

template<typename T>
inline string to_string(const af::array a, string_view fmt = "%s") {
    auto b = af::flat(a);
    b.eval();
    string s = "[";
    for (auto i : range(b.elements())) {
        if (i > 0) s += ' ';
        // TODO command in format() to auto add separator (if not already present and if not start of string)
        format_s(s, fmt, b(i).scalar<T>());
    }
    s += ']';
    return s;
}

inline auto back(af::dim4 a) { return a[a.ndims() - 1]; }
inline auto pop_front(af::dim4 a) { return af::dim4(a[1], a[2], a[3], 1); }
