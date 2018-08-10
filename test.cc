#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "callstack.h"
#include "timestamp.h"

void sigsegv_handler(int sig) {
	fprintf(stderr, "Error: signal %d:\n", sig);
	Callstack stack;
	std::string s;
	stack.write(s, {"sigsegv_handler(int)", "_sigtramp"});
	fputs(s.c_str(), stderr);
	exit(1);
}

int main(int argc, char** argv) {
	void sigsegv_handler(int sig);
	signal(SIGSEGV, sigsegv_handler);
	Timestamp::init();
	return Catch::Session().run(argc, argv);
}
