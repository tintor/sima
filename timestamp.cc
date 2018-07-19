#include "timestamp.h"
#include "util.h"
#include <chrono>
#include <unistd.h>
#include <iostream>

void Timestamp::init() {
	auto p = std::chrono::high_resolution_clock::period();
	std::cout << "HRC " << double(p.num) / p.den << std::endl;

	auto ax = std::chrono::high_resolution_clock::now();
	long bx = __builtin_readcyclecounter();
	usleep(100000);
	auto ay = std::chrono::high_resolution_clock::now();
	long by = __builtin_readcyclecounter();
	usleep(100000);
	auto az = std::chrono::high_resolution_clock::now();
	long bz = __builtin_readcyclecounter();

	usleep(100000);

        auto cx = std::chrono::high_resolution_clock::now();
        long dx = __builtin_readcyclecounter();
    	usleep(100000);
        auto cy = std::chrono::high_resolution_clock::now();
        long dy = __builtin_readcyclecounter();
    	usleep(100000);
    	auto cz = std::chrono::high_resolution_clock::now();
        long dz = __builtin_readcyclecounter();

        std::chrono::duration<double, std::milli> cax = cx - ax;
        std::chrono::duration<double, std::milli> cay = cy - ay;
        std::chrono::duration<double, std::milli> caz = cz - az;

        long dbx = dx - bx;
        long dby = dy - by;
        long dbz = dz - bz;

    	s_milisec_per_tick = median(cax.count() / dbx, cay.count() / dby, caz.count() / dbz) * 1000;
        std::cout << "mpt " << s_milisec_per_tick << std::endl;
}

double Timestamp::s_milisec_per_tick = 0;
