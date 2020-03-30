#include <core/callstack.h>
#include <core/string_util.h>
#include <core/auto.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <signal.h>

using namespace std::literals;

Callstack::Callstack() {
	_size = backtrace(_stack.data(), _stack.size());
}

std::unique_ptr<char, void(*)(char*)> Callstack::demangle(const char* symbol) {
	int status;
	char* e = abi::__cxa_demangle(symbol, nullptr, nullptr, &status);
	return {(status == 0) ? e : nullptr, [](char* p) { free(p); }};
}

static bool matches(const std::initializer_list<string_view>& prefixes, string_view s) {
	for (string_view m : prefixes)
		if (s.size() >= m.size() && s.substr(0, m.size()) == m)
			return true;
	return false;
}

void Callstack::write(string& out, const std::initializer_list<string_view>& exclude) const {
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
		if (status != 0 || (buffer != "Callstack::Callstack()"sv && !matches(exclude, buffer))) {
			for (int j = 0; j < sp.size(); j++) {
				out += (j == 3 && status == 0) ? string_view(buffer) : sp[j];
				out += ' ';
			}
			out += '\n';
		}
	}
}

static void sigsegv_handler(int sig) {
	fprintf(stderr, "Error: signal %d:\n", sig);
	Callstack stack;
	string s;
	stack.write(s, {"sigsegv_handler(int)", "_sigtramp"});
	fputs(s.c_str(), stderr);
	exit(1);
}

void InitSegvHandler() {
	signal(SIGSEGV, sigsegv_handler);
}
