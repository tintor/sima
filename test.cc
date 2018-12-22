#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <core/callstack.h>
#include <core/timestamp.h>

int main(int argc, char** argv) {
	InitSegvHandler();
	Timestamp::init();
	return Catch::Session().run(argc, argv);
}
