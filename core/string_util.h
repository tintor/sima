#pragma once
#include <core/std.h>
#include <core/util.h>
#include <strstream>

template<typename Number>
Number parse(string_view s) {
    std::istrstream is(s.data(), s.size());
	Number out;
    is >> out;
	return out;
}

template<typename Number>
Number parse(pair<const char*, const char*> s) {
    std::istrstream is(s.first, s.second - s.first);
	Number out;
    is >> out;
	return out;
}

inline bool search(string_view s, const regex& re) {
	return std::regex_search(s.begin(), s.end(), re);
}

inline bool search(string_view s, const regex& re, std::cmatch& match) {
	return std::regex_search(s.begin(), s.end(), /*out*/match, re);
}

inline bool match(string_view s, const regex& re) {
	return std::regex_match(s.begin(), s.end(), re);
}

inline bool match(string_view s, const regex& re, std::cmatch& match) {
	return std::regex_match(s.begin(), s.end(), /*out*/match, re);
}

inline vector<string_view> split(string_view s, cspan<char> delim) {
	vector<string_view> out;
	const char* b = s.begin();
	int c = 0;
	for (const char* i = s.begin(); i != s.end(); i++) {
		if (contains(delim, *i)) {
			if (c > 0)
				out.push_back(string_view(b, c));
			b = i + 1;
			c = 0;
		} else {
			c += 1;
		}
	}
	if (c > 0)
		out.push_back(string_view(b, c));
	return out;
}

inline vector<string_view> split(string_view s, char delim = ' ') {
	cspan<char> d = { delim };
	return split(s, d);
}

inline bool is_digit(char c) {
	return '0' <= c && c <= '9';
}

// TODO make it work for large numbers
inline bool natural_less(string_view a, string_view b) {
	auto ai = a.begin(), bi = b.begin();
	while (ai != a.end() && bi != b.end()) {
		bool an = false, bn = false;
		size_t av, bv;

		if (is_digit(*ai)) {
			an = true;
			av = 0;
			while (ai != a.end() && is_digit(*ai)) {
				av = av * 10 + (*ai - '0');
				ai++;
			}
		} else
			av = *ai++;

		if (is_digit(*bi)) {
			bn = true;
			bv = 0;
			while (bi != b.end() && is_digit(*bi)) {
				bv = bv * 10 + (*bi - '0');
				bi++;
			}
		} else
			bv = *bi++;

		if (an && !bn)
			av = '0';
		if (!an && bn)
			bv = '0';
		if (av != bv)
			return av < bv;
	}
	return ai == a.end() && bi != b.end();
}

inline string cat(string_view a, string_view b) {
	string s;
	s.resize(a.size() + b.size());
	copy(a.begin(), a.end(), s.begin());
	copy(b.begin(), b.end(), s.begin() + a.size());
	return s;
}

inline string cat(string_view a, const string& b) {
	string s;
	s.resize(a.size() + b.size());
	copy(a.begin(), a.end(), s.begin());
	copy(b.begin(), b.end(), s.begin() + a.size());
	return s;
}
