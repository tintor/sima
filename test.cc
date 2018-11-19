#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "callstack.h"
#include "timestamp.h"

int main(int argc, char** argv) {
	InitSegvHandler();
	Timestamp::init();
	return Catch::Session().run(argc, argv);
}
