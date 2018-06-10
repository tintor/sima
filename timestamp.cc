#include "timestamp.hh"
#include "util.hh"
#include <chrono>
#include <unistd.h>
#include <iostream>

void Timestamp::init() {
	auto p = std::chrono::high_resolution_clock::period();
	std::cout << "HRC " << double(p.num) / p.den << std::endl;

    auto ax = std::chrono::high_resolution_clock::now();
    int64_t bx = rdtsc();
    	usleep(100000);
        auto ay = std::chrono::high_resolution_clock::now();
        int64_t by = rdtsc();
    	usleep(100000);
    	auto az = std::chrono::high_resolution_clock::now();
        int64_t bz = rdtsc();

        usleep(100000);

        auto cx = std::chrono::high_resolution_clock::now();
        int64_t dx = rdtsc();
    	usleep(100000);
        auto cy = std::chrono::high_resolution_clock::now();
        int64_t dy = rdtsc();
    	usleep(100000);
    	auto cz = std::chrono::high_resolution_clock::now();
        int64_t dz = rdtsc();

        std::chrono::duration<double, std::milli> cax = cx - ax;
        std::chrono::duration<double, std::milli> cay = cy - ay;
        std::chrono::duration<double, std::milli> caz = cz - az;

        int64_t dbx = dx - bx;
        int64_t dby = dy - by;
        int64_t dbz = dz - bz;

    	s_milisec_per_tick = median(cax.count() / dbx, cay.count() / dby, caz.count() / dbz) * 1000;
        std::cout << "mpt " << s_milisec_per_tick << std::endl;
}

double Timestamp::s_milisec_per_tick = 0;
