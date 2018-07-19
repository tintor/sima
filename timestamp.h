#pragma once

#include <chrono>

/*inline int64_t rdtsc() {
	uint32_t lo, hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return (static_cast<uint64_t>(hi) << 32) | lo;
}*/

struct Timestamp {
	Timestamp() : m_ticks(__builtin_readcyclecounter()) { }

	long elapsed(Timestamp a = Timestamp()) const { return a.m_ticks - m_ticks; }

	double elapsed_ms(Timestamp a = Timestamp()) const { return elapsed(a) * s_milisec_per_tick; }

	static void init();

private:
	long m_ticks;
	static double s_milisec_per_tick;
};
