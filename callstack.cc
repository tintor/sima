#include "callstack.h"
#include "util.h"
#include "auto.h"
#include <execinfo.h>
#include <cxxabi.h>

Callstack::Callstack() {
	_size = backtrace(_stack.data(), _stack.size());
}

std::unique_ptr<char, void(*)(char*)> Callstack::demangle(const char* symbol) {
	int status;
	char* e = abi::__cxa_demangle(symbol, nullptr, nullptr, &status);
	return {(status == 0) ? e : nullptr, [](char* p) { free(p); }};
}

void Callstack::write(std::string& out) const {
	char** strs = backtrace_symbols(_stack.data(), _size);
	ON_SCOPE_EXIT(free(strs));

	size_t length = 256;
	char* buffer = (char*)malloc(length);
	ON_SCOPE_EXIT(free(buffer));

	for (int i = 0; i < _size; ++i) {
		auto sp = split(strs[i]);
		char* symbol = const_cast<char*>(sp[3].data());
		symbol[sp[3].size()] = '\0';
		int status;
		char* new_buffer = abi::__cxa_demangle(symbol, buffer, &length, &status);
		if (new_buffer)
			buffer = new_buffer;
		if (status != 0 || buffer != std::string_view("Callstack::Callstack()")) {
			for (int j = 0; j < sp.size(); j++) {
				out += (j == 3 && status == 0) ? std::string_view(buffer) : sp[j];
				out += ' ';
			}
			out += '\n';
		}
	}
}
