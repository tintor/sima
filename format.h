#pragma once
// Type-safe and extensible C++17 string formatting library:
// std::string s;
// format(s, "Hello {} world!"sv, 20);
//
// It also adds operator<< for std::string similar to ostream (but faster as
// ostream is heavy).
#include <array>
#include <string>
#include <string_view>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>

template <typename... Args>
void format_s(std::string& s, std::string_view fmt, Args... args);

template<typename T>
void format_e(std::string& s, std::string_view spec, const std::vector<T>& m) {
	s += '[';
	if (m.size() > 0) {
		format_s(s, "%s", m[0]);
		for (size_t i = 1; i < m.size(); i++)
			format_s(s, ", %s", m[i]);
	}
	s += ']';
}

template<typename T, size_t S>
void format_e(std::string& s, std::string_view spec, const std::array<T, S>& m) {
	s += '[';
	if (m.size() > 0) {
		format_s(s, "%s", m[0]);
		for (size_t i = 1; i < m.size(); i++)
			format_s(s, ", %s", m[i]);
	}
	s += ']';
}

template<typename Map>
void format_map(std::string& s, std::string_view spec, const Map& m) {
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

#define FORMAT_MAP(MAP) template<typename Key, typename Value, typename Compare> inline void format_e(std::string& s, std::string_view spec, const MAP<Key, Value, Compare>& m) { format_map(s, spec, m); }

FORMAT_MAP(std::map);
FORMAT_MAP(std::multimap);
FORMAT_MAP(std::unordered_map);
FORMAT_MAP(std::unordered_multimap);

template<typename Set>
void format_set(std::string& s, std::string_view spec, const Set& m) {
	s += '{';
	bool first = true;
	if (m.size() > 0) {
		for (const auto& e : m) {
			if (!first)
				s += ", ";
			format_s(s, "%s", e);
			first = false;
		}
	}
	s += '}';
}

#define FORMAT_SET(SET) template<typename Key, typename Compare> inline void format_e(std::string& s, std::string_view spec, const SET<Key, Compare>& m) { format_set(s, spec, m); }

FORMAT_SET(std::set);
FORMAT_SET(std::multiset);
FORMAT_SET(std::unordered_set);
FORMAT_SET(std::unordered_multiset);

template<typename A, typename B>
inline void format_e(std::string& s, std::string_view spec, const std::pair<A, B>& p) {
	using namespace std::literals;
	format_s(s, "(%s %s)"sv, p.first, p.second);
}

inline void format_e(std::string& s, std::string_view spec, const std::string& v) { s += v; }

inline void format_e(std::string& s, std::string_view spec, std::string_view v) { s += v; }

inline void format_e(std::string& s, std::string_view spec, const char* v) { s += v ? v : "null"; }

inline void format_e(std::string& s, std::string_view spec, char v) { s += v; }

inline void format_e(std::string& s, std::string_view spec, bool v) {
	using namespace std::literals;
	s += v ? "true"sv : "false"sv;
}

template <typename Integer>
void format_int(std::string& s, std::string_view spec, Integer v) {
	if (v == 0) {
		s += '0';
		return;
	}
	// TODO hexadecimal if spec == '%x'
	std::array<char, /*sign*/ 1 + sizeof(v) * 5 / 2> buffer;
	char* p = buffer.end();
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
	s += std::string_view(p, buffer.end() - p);
}

#define FORMAT_INT(T)                                                          \
	inline void format_e(std::string& s, std::string_view spec, T v) {         \
		format_int(s, spec, v);                                                \
	}

FORMAT_INT(short)
FORMAT_INT(unsigned short)
FORMAT_INT(int)
FORMAT_INT(unsigned int)
FORMAT_INT(long)
FORMAT_INT(unsigned long)
FORMAT_INT(__int128)
FORMAT_INT(__uint128_t)

using cent = __int128;

inline void format_e(std::string& s, std::string_view spec, float v) {
	// TODO pass spec to snprintf
	std::array<char, 30> buffer;
	int length = std::snprintf(buffer.data(), buffer.size(), "%f", v);
	if (length <= 0)
		throw std::runtime_error("format: snprintf failed");
	s += std::string_view(buffer.data(), length);
}

inline void format_e(std::string& s, std::string_view spec, double v) {
	std::array<char, 30> buffer;
	int length = std::snprintf(buffer.data(), buffer.size(), "%lf", v);
	if (length <= 0)
		throw std::runtime_error("format: snprintf failed");
	s += std::string_view(buffer.data(), length);
}

inline void format_e(std::string& s, std::string_view spec, long double v) {
	std::array<char, 60> buffer;
	int length = std::snprintf(buffer.data(), buffer.size(), "%Lf", v);
	if (length <= 0)
		throw std::runtime_error("format: snprintf failed");
	s += std::string_view(buffer.data(), length);
}

#define FORMAT_VEC(T, N) \
using T##N = T __attribute__((ext_vector_type(N))); \
inline void format_e(std::string& s, std::string_view spec, T##N v) { \
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

/*inline void format_e(std::string& s, std::string_view spec, int __attribute__((ext_vector_type(2))) v) {
	s += '[';
	format_e(s, "", v.x);
	s += ' ';
	format_e(s, "", v.y);
	s += ']';
}

inline void format_e(std::string& s, std::string_view spec, float __attribute__((ext_vector_type(2))) v) {
	s += '[';
	format_e(s, "", v.x);
	s += ' ';
	format_e(s, "", v.y);
	s += ']';
}

template<typename T>
inline void format_e(std::string& s, std::string_view spec, T __attribute__((ext_vector_type(3))) v) {
	s += '[';
	format_e(s, "", v.x);
	s += ' ';
	format_e(s, "", v.y);
	s += ' ';
	format_e(s, "", v.z);
	s += ']';
}

template<typename T>
inline void format_e(std::string& s, std::string_view spec, T __attribute__((ext_vector_type(4))) v) {
	s += '[';
	format_e(s, "", v.x);
	s += ' ';
	format_e(s, "", v.y);
	s += ' ';
	format_e(s, "", v.z);
	s += ' ';
	format_e(s, "", v.w);
	s += ']';
}

template<typename T>
inline void format_e(std::string& s, std::string_view spec, T __attribute__((ext_vector_type(8))) v) {
	s += '[';
	format_e(s, "", v.x);
	for (int i = 1; i < 8; i++) {
		s += ' ';
		format_e(s, "", v[i]);
	}
	s += ']';
}*/

inline void format_a(std::string& s, int skip, std::string_view spec) {
	throw std::runtime_error("format: not enough arguments");
}

template <typename T, typename... Args>
void format_a(std::string& s, int skip, std::string_view spec, const T& arg,
			  Args... args) {
	if (skip == 0) {
		format_e(s, spec, arg);
	} else {
		format_a(s, skip - 1, spec, args...);
	}
}

// TODO check if string_view::find() is faster
// TODO send multiple chars at once to s instead of one by one
// TODO support positional args
template <typename... Args>
void format_s(std::string& s, std::string_view fmt, Args... args) {
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
			if (p == e)
				throw std::runtime_error("format: unterminated %");
			c = *p++;
			if (c == '%' && p - q == 1) {
				s += '%';
				break;
			}
			if (c == 's' || c == 'f' || c == 'd' || c == 'g' || c == 'x' || c == 'X' || c == 'p') {
				format_a(s, a++, std::string_view(q, p - q - 1), args...);
				break;
			}
		}
	}
}

template <typename T> std::string& operator<<(std::string& s, T v) {
	using namespace std::literals;
	format_e(s, ""sv, v);
	return s;
}

template <typename... Args>
std::string format(std::string_view fmt, Args... args) {
	std::string s;
	format_s(s, fmt, args...);
	return s;
}

template <typename... Args>
void print(std::string_view fmt, Args... args) {
	std::cout << format(fmt, args...);
}

#ifdef NDEBUG
template <typename... Args>
void dprint(std::string_view fmt, Args... args) {
}
#else
template <typename... Args>
void dprint(std::string_view fmt, Args... args) {
	std::cout << format(fmt, args...);
}
#endif
