#pragma once
#include "std.h"
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

inline vector<string_view> split(string_view s, char delim = ' ') {
	vector<string_view> out;
	const char* b = s.begin();
	int c = 0;
	for (const char* i = s.begin(); i != s.end(); i++) {
		if (*i == delim) {
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
