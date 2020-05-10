#pragma once
#include <core/std.h>

template <typename T>
class matrix {
   private:
    uint _dim_a = 0;
    uint _dim_b = 0;
    vector<T> _data;

   public:
    void resize(uint a, uint b) {
        _dim_a = a;
        _dim_b = b;
        _data.resize(a * b);
    }

    void fill(T value) {
        for (T& e : _data) e = value;
    }

    const T& operator()(uint a, uint b) const {
        if (a >= _dim_a || b >= _dim_b) THROW(invalid_argument);
        return _data[b * _dim_a + a];
    }

    T& operator()(uint a, uint b) {
        if (a >= _dim_a || b >= _dim_b) THROW(invalid_argument);
        return _data[b * _dim_a + a];
    }
};
