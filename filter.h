#pragma once
#include "std.h"

template<typename T>
class LowPassFilter {
public:
	LowPassFilter(T param, T init) : m_param(param), m_accum(init) { }

	T tick(T a) {
		m_accum = a * m_param + m_accum * (1 - m_param);
		return m_accum;
	}

private:
	const T m_param;
	T m_accum;
};

template<typename T>
class LagFilter {
public:
	LagFilter(uint lag, T init) {
		m_buffer.resize(lag, init);
	}

	T tick(T a) {
		if (m_buffer.size() == 0)
			return a;
		if (m_buffer.size() == 1) {
			T e = m_buffer[0];
			m_buffer[0] = a;
			return e;
		}
		uint i = m_pos + 1;
		if (i >= m_buffer.size())
			i -= m_buffer.size();
		T e = m_buffer[i];
		m_buffer[m_pos] = a;
		m_pos += 1;
		if (m_pos == m_buffer.size())
			m_pos = 0;
		return e;
	}

private:
	uint m_pos = 0;
	vector<T> m_buffer;
};
