#pragma once

#include <strstream>
#include <algorithm>
#include <regex>

template<typename T>
T min(std::initializer_list<T> args) {
	T result = std::numeric_limits<T>::max();
	for (auto e : args)
		result = std::min(e, result);
	return result;
}

template<typename T>
T max(std::initializer_list<T> args) {
	T result = std::numeric_limits<T>::min();
	for (auto e : args)
		result = std::max(e, result);
	return result;
}

template<typename T>
T median(T x, T y, T z) {
	if (x <= y) {
		if (z <= x)
			return x;
		if (y <= z)
			return y;
		return z;
	}
	// y < x
	if (z <= y)
		return y;
	if (x <= z)
		return x;
	return z;
}

template<typename Container>
void release(Container& container) {
	Container empty;
	std::swap(empty, container);
}

template<typename Container>
void sort(Container& container) {
	std::sort(container.begin(), container.end());
}

template<typename Container, typename Func>
void sort(Container& container, const Func& func) {
	std::sort(container.begin(), container.end(), func);
}

// from C++17, but not in OSX clang yet
template<typename Number>
void from_chars(const char* begin, const char* end, Number& out) {
    std::istrstream is(begin, end - begin);
    is >> out;
}

template<typename Number>
void from_chars(std::string_view s, Number& out) {
    std::istrstream is(s.data(), s.size());
    is >> out;
}

template<typename Number>
void from_chars(std::pair<const char*, const char*> s, Number& out) {
    std::istrstream is(s.first, s.second - s.first);
    is >> out;
}

template<typename Number>
Number parse(std::string_view s) {
    std::istrstream is(s.data(), s.size());
	Number out;
    is >> out;
	return out;
}

template<typename Number>
Number parse(std::pair<const char*, const char*> s) {
    std::istrstream is(s.first, s.second - s.first);
	Number out;
    is >> out;
	return out;
}

inline bool search(std::string_view s, const std::regex& re) {
	return std::regex_search(s.begin(), s.end(), re);
}		

inline bool search(std::string_view s, const std::regex& re, std::cmatch& match) {
	return std::regex_search(s.begin(), s.end(), /*out*/match, re);
}		

inline bool match(std::string_view s, const std::regex& re) {
	return std::regex_match(s.begin(), s.end(), re);
}		

inline bool match(std::string_view s, const std::regex& re, std::cmatch& match) {
	return std::regex_match(s.begin(), s.end(), /*out*/match, re);
}		

class error : public std::runtime_error {
public:
	template <typename... Args>
	error(std::string_view fmt, Args... args) : std::runtime_error(format(fmt, args...)) {
	}
};

template<typename Vector>
void remove_dups(Vector& v) {
	std::sort(v.begin(), v.end());
	v.erase(std::unique(v.begin(), v.end()), v.end());
}

inline std::vector<std::string_view> split(std::string_view s, char delim = ' ') {
	std::vector<std::string_view> out;
	const char* b = s.begin();
	int c = 0;
	for (const char* i = s.begin(); i != s.end(); i++) {
		if (*i == delim) {
			if (c > 0)
				out.push_back(std::string_view(b, c));
			b = i + 1;
			c = 0;
		} else {
			c += 1;
		}
	}
	if (c > 0)
		out.push_back(std::string_view(b, c));
	return out;
}

