#define CATCH_CONFIG_RUNNER
#include <core/callstack.h>
#include <core/timestamp.h>

#include <catch.hpp>

int main(int argc, char** argv) {
    InitSegvHandler();
    Timestamp::init();
    return Catch::Session().run(argc, argv);
}
