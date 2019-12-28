#pragma once
// Type-safe and extensible C++17 string formatting library:
// string s;
// format(s, "Hello {} world!"sv, 20);
//
// It also adds operator<< for string similar to ostream (but faster as
// ostream is heavy).
#include <core/int.h>
#include <core/std.h>
#include <core/span.h>
#include <core/interval.h>
#include <complex>

template<typename... Args>
void format_s(string& s, string_view fmt, const Args& ... args);

template<typename T>
void format_span(string& s, string_view spec, span<T> m) {
	s += '[';
	if (m.size() > 0) {
		format_s(s, spec, m[0]);
		for (size_t i = 1; i < m.size(); i++) {
			s += ", ";
			format_s(s, spec, m[i]);
		}
	}
	s += ']';
}

template<typename T>
void format_e(string& s, string_view spec, span<T> m) {
	format_span(s, spec, m);
}

template<typename T, typename A>
void format_e(string& s, string_view spec, const vector<T, A>& m) {
	format_span(s, spec, cspan<T>(m.data(), m.size()));
}

template<size_t S>
void format_e(string& s, string_view spec, const array<char, S>& m) {
	s.append(m.data(), strnlen(m.data(), m.size()));
}

template<typename T, size_t S>
void format_e(string& s, string_view spec, const array<T, S>& m) {
	format_span(s, spec, span<const T>(m));
}

template<typename Map>
void format_map(string& s, string_view spec, const Map& m) {
	s += '{';
	bool first = true;
	if (m.size() > 0) {
		for (const auto& [k, v] : m) {
			if (!first)
				s += ", ";
			format_s(s, "%s: %s", k, v);
			first = false;
		}
	}
	s += '}';
}

#define FORMAT_MAP(MAP) template<typename Key, typename Value, typename Compare> \
	inline void format_e(string& s, string_view spec, const MAP<Key, Value, Compare>& m) { format_map(s, spec, m); }

FORMAT_MAP(map);
FORMAT_MAP(multimap);
FORMAT_MAP(unordered_map);
FORMAT_MAP(unordered_multimap);

template<typename Set>
void format_set(string& s, string_view spec, const Set& m) {
	s += '{';
	bool first = true;
	if (m.size() > 0) {
		for (const auto& e : m) {
			if (!first)
				s += ", ";
			format_s(s, spec, e);
			first = false;
		}
	}
	s += '}';
}

#define FORMAT_SET(SET) template<typename Key, typename Compare> \
	inline void format_e(string& s, string_view spec, const SET<Key, Compare>& m) { format_set(s, spec, m); }

FORMAT_SET(set);
FORMAT_SET(multiset);
FORMAT_SET(unordered_set);
FORMAT_SET(unordered_multiset);

template<typename A, typename B>
inline void format_e(string& s, string_view spec, const pair<A, B>& p) {
	format_s(s, "(%s %s)", p.first, p.second);
}

template<typename A>
inline void format_e(string& s, string_view spec, const optional<A>& p) {
	if (p.has_value())
		format_s(s, "%s", *p);
	else
		s += "null";
}

inline void format_e(string& s, string_view spec, const string& v) { s += v; }
inline void format_e(string& s, string_view spec, string_view v) { s += v; }
inline void format_e(string& s, string_view spec, const char* v) { s += v ? v : "null"; }
inline void format_e(string& s, string_view spec, char v) { s += v; }
inline void format_e(string& s, string_view spec, bool v) { s += v ? "true"sv : "false"sv; }

template<typename Integer>
void format_int(string& s, Integer v, int zero_pad, bool hexadecimal) {
	if (v == 0) {
		zero_pad -= 1;
		for (int i = 0; i < zero_pad; i++)
			s += '0';

		s += '0';
		return;
	}

	array<char, /*sign*/ 1 + sizeof(v) * 5 / 2> buffer;
	char* p = buffer.end();
	const char* digits = "0123456789ABCDEF";

	if (hexadecimal) {
		if (v < 0) {
			while (v != 0) {
				assert(buffer.begin() < p);
				*--p = digits[-(v % 16)];
				v /= 16;
			}
			assert(buffer.begin() < p);
			*--p = '-';
		} else {
			while (v != 0) {
				assert(buffer.begin() < p);
				*--p = digits[v % 16];
				v /= 16;
			}
		}
	} else {
		if (v < 0) {
			while (v != 0) {
				assert(buffer.begin() < p);
				*--p = '0' - (v % 10);
				v /= 10;
			}
			assert(buffer.begin() < p);
			*--p = '-';
		} else {
			while (v != 0) {
				assert(buffer.begin() < p);
				*--p = '0' + (v % 10);
				v /= 10;
			}
		}
	}
	zero_pad -= buffer.end() - p;
	for (int i = 0; i < zero_pad; i++)
		s += '0';
	s += string_view(p, buffer.end() - p);
}

template<typename T>
inline T abs(T a, T b) {
	return (a > b) ? a - b : b - a;
}

inline bool try_short(string& s, ucent v, ulong m, char c) {
	constexpr int E = 100;
	if (v < m - m / E)
		return false;

	auto a = (v + m / 2) / m;
	if (abs(a * m, v) * E > v)
		return false;
	format_int(s, a, 0, false);
	s += c;
	return true;
}

template<typename Integer>
void format_int(string& s, string_view spec, Integer v) {
	if (spec == "%h"sv && v != -v) {
		// print int with order of magnitude, ie 123m
		if (v < 0) {
			s += '-';
			v = -v;
		}

		if (try_short(s, v, 1000000000000lu, 't'))
			return;
		if (try_short(s, v, 1000000000lu, 'g'))
			return;
		if (try_short(s, v, 1000000, 'm'))
			return;
		if (try_short(s, v, 1000, 'k'))
			return;
		format_int(s, v, 0, false);
		return;
	}

	if (spec == "%t"sv && v != -v) {
		// interpret v as time in seconds: and print as [[hh:]mm:]ss
		if (v < 0) {
			s += '-';
			v = -v;
		}

		uint sec = v % 60;
		v /= 60;
		uint min = v % 60;
		auto hour = v / 60;

		if (hour != 0) {
			format_int(s, hour, 0, false);
			s += ':';
			format_int(s, min, 2, false);
			s += ':';
			format_int(s, sec, 2, false);
		} else if (min != 0) {
			format_int(s, min, 0, false);
			s += ':';
			format_int(s, sec, 2, false);
		} else
			format_int(s, sec, 0, false);
		return;
	}

	int zero_pad = 0;
	if (spec.size() >= 3 && spec[0] == '%' && spec[1] == '0' && '1' <= spec[2] && spec[2] <= '9')
		zero_pad = spec[2] - '0';

	bool hexadecimal = spec.size() >= 2 && spec[0] == '%' && spec.back() == 'x';

	format_int(s, v, zero_pad, hexadecimal);
}

#define FORMAT_INT(T) inline void format_e(string& s, string_view spec, const T& v) { format_int(s, spec, v); }

FORMAT_INT(short)
FORMAT_INT(ushort)
FORMAT_INT(int)
FORMAT_INT(uint)
FORMAT_INT(long)
FORMAT_INT(ulong)
FORMAT_INT(cent)
FORMAT_INT(ucent)

template<typename T>
inline void format_e(string& s, string_view spec, const T* v) { format_int(s, "%x", reinterpret_cast<size_t>(v)); }

// TODO precompute
inline double POW10(int a) { return std::pow(10, a); }

inline void format_fp(string& s, double fvalue, int min, int max, bool plus, bool minus, bool zero, bool space) {
        int signvalue = 0;
        char iconvert[311];
        char fconvert[311];
        int iplace = 0;
        int fplace = 0;
        int padlen = 0; /* amount to pad */
        int zpadlen = 0;
        int index;
        double intpart;
        double fracpart;
        double temp;

        if (max < 0)
                max = 6;

        double ufvalue = std::abs(fvalue);

        if (fvalue < 0) {
                signvalue = '-';
        } else {
                if (plus) {
                        signvalue = '+';
                } else {
                        if (space)
                                signvalue = ' ';
                }
        }

        /*
         * Sorry, we only support 16 digits past the decimal because of our
         * conversion method
         */
        if (max > 16)
                max = 16;

        /* We "cheat" by converting the fractional part to integer by
         * multiplying by a factor of 10
         */

        temp = ufvalue;
		std::modf(temp, &intpart);

        fracpart = std::round((POW10(max)) * (ufvalue - intpart));

        if (fracpart >= POW10(max)) {
                intpart++;
                fracpart -= POW10(max);
        }

        // Convert integer part
        do {
                temp = intpart;
                std::modf(intpart * 0.1, &intpart);
                temp = temp * 0.1;
                index = (int) ((temp - intpart + 0.05) * 10.0);
                iconvert[iplace++] = '0' + index;
        } while (intpart && (iplace < 311));
        if (iplace == 311) iplace--;
        iconvert[iplace] = 0;

        /* Convert fractional part */
        if (fracpart) {
                do {
                        temp = fracpart;
						std::modf(fracpart * 0.1, &fracpart);
                        temp = temp * 0.1;
                        index = (int) ((temp - fracpart + 0.05) * 10.0);
                        fconvert[fplace++] = '0' + index;
                } while(fracpart && (fplace < 311));
                if (fplace == 311) fplace--;
        }
        fconvert[fplace] = 0;

        /* -1 for decimal point, another -1 if we are printing a sign */
        padlen = min - iplace - max - 1 - (signvalue ? 1 : 0);
        zpadlen = max - fplace;
        if (zpadlen < 0)
				zpadlen = 0;
        if (padlen < 0)
                padlen = 0;
        if (minus)
                padlen = -padlen; /* Left Justifty */
        if (zero && (padlen > 0)) {
                if (signvalue) {
						s += signvalue;
                        --padlen;
                        signvalue = 0;
                }
                while (padlen > 0) {
                        s += '0';
                        --padlen;
                }
        }
        while (padlen > 0) {
				s += ' ';
                --padlen;
        }
        if (signvalue)
                s += signvalue;
        while (iplace > 0)
                s += iconvert[--iplace];

        // Decimal point.
        if (max > 0) {
				s += '.';
                while (fplace > 0)
                        s += fconvert[--fplace];
        }
        while (zpadlen > 0) {
				s += '0';
                --zpadlen;
        }
	while (padlen < 0) {
		s += ' ';
		++padlen;
	}
}

inline bool isIntSpec(string_view spec) {
	return spec.size() >= 2 && spec[0] == '%' && (spec.back() == 'd' || spec.back() == 'b' || spec.back() == 'x' || spec.back() == 'b' || spec.back() == 't');
}

template<typename T>
inline void formatFloatAsInt(string& s, string_view spec, T v) {
	if (v < 0)
			format_int(s, spec, long(std::round(v)));
		else
			format_int(s, spec, ulong(std::round(v)));
}

inline void format_e(string& s, string_view spec, float v) {
	if (isIntSpec(spec)) {
		formatFloatAsInt(s, spec, v);
		return;
	}
	array<char, 30> buffer;
	// TODO don't use string
	// TODO reset spec to f if lf
	string _spec(spec);
	int length = std::snprintf(buffer.data(), buffer.size(), spec == "%s" ? "%g" : _spec.c_str(), v);
	if (length <= 0)
		throw std::runtime_error("format: snprintf failed");
	s += string_view(buffer.data(), length);
}

inline void format_e(string& s, string_view spec, double v) {
	if (isIntSpec(spec)) {
		formatFloatAsInt(s, spec, v);
		return;
	}
	array<char, 30> buffer;
	string _spec(spec);
	int length = std::snprintf(buffer.data(), buffer.size(), spec == "%s" ? "%g" : _spec.c_str(), v);
	if (length <= 0)
		throw std::runtime_error("format: snprintf failed");
	s += string_view(buffer.data(), length);
}

inline void format_e(string& s, string_view spec, long double v) {
	if (isIntSpec(spec)) {
		formatFloatAsInt(s, spec, v);
		return;
	}
	array<char, 60> buffer;
	string _spec(spec);
	int length = std::snprintf(buffer.data(), buffer.size(), spec == "%s" ? "%Lg" : _spec.c_str(), v);
	if (length <= 0)
		throw std::runtime_error("format: snprintf failed");
	s += string_view(buffer.data(), length);
}

template<typename T>
void format_e(string& s, string_view spec, complex<T> v) {
	double re = v.real(), im = v.imag();
	if (im > 0) {
		if (re != 0)
			format_e(s, spec, re);
		s += '+';
		format_e(s, spec, im);
		s += 'i';
		return;
	}
	if (im < 0) {
		if (re != 0)
			format_e(s, spec, re);
		format_e(s, spec, im);
		s += 'i';
		return;
	}
	format_e(s, spec, re);
}

#define FORMAT_VEC(T, N) \
using T##N = T __attribute__((ext_vector_type(N))); \
inline void format_e(string& s, string_view spec, T##N v) { \
	format_e(s, spec, v.x); \
	for (int i = 1; i < N; i++) { \
		s += ' '; \
		format_e(s, spec, v[i]); \
	} \
}

#define FORMAT_VECN(T) FORMAT_VEC(T, 2); FORMAT_VEC(T, 3); FORMAT_VEC(T, 4); FORMAT_VEC(T, 8)

FORMAT_VECN(int);
FORMAT_VECN(long);
FORMAT_VECN(cent);
FORMAT_VECN(float);
FORMAT_VECN(double);

template<typename T>
void format_e(string& s, string_view spec, interval<T> m) {
	s += '[';
	format_e(s, spec, m.min);
	s += ' ';
	format_e(s, spec, m.max);
	s += ']';
}

inline void format_a(string& s, int skip, string_view spec) {
	throw std::runtime_error("format: not enough arguments");
}

template<typename T, typename... Args>
void format_a(string& s, int skip, string_view spec, const T& arg,
			  const Args& ... args) {
	if (skip == 0) {
		format_e(s, spec, arg);
	} else {
		format_a(s, skip - 1, spec, args...);
	}
}

template<typename T>
void format_e(string& s, string_view spec, const atomic<T>& m) {
	format_e(s, spec, m.load());
}

// TODO check if string_view::find() is faster
// TODO send multiple chars at once to s instead of one by one
// TODO support positional args
template<typename... Args>
void format_s(string& s, string_view fmt, const Args& ... args) {
	auto p = fmt.begin(), e = fmt.end();
	int a = 0;
	while (p != e) {
		char c = *p++;
		if (c != '%') {
			s += c;
			continue;
		}
		auto q = p;
		while (true) {
			if (p == e) {
				string m = "format: unterminated % in [";
				m += fmt;
				m += ']';
				throw std::runtime_error(m);
			}
			c = *p++;
			if (c == '%' && p - q == 1) {
				s += '%';
				break;
			}
			if (c == 's' || c == 'f' || c == 'd' || c == 'g' || c == 'x' || c == 'X' || c == 'p' || c == 't' || c == 'h') {
				format_a(s, a++, string_view(q - 1, p - q + 1), args...);
				break;
			}
		}
	}
}

template<typename T> string& operator<<(string& s, T v) {
	format_e(s, "", v);
	return s;
}

template<typename... Args>
[[nodiscard]] string format(string_view fmt, const Args& ... args) {
	string s;
	format_s(s, fmt, args...);
	return s;
}

template<typename... Args>
void print(string_view fmt, const Args& ... args) {
	cout << format(fmt, args...);
}

template<typename... Args>
void dprint(string_view fmt, const Args& ... args) {
	if (debug)
		cout << format(fmt, args...);
}
