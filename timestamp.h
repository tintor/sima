#pragma once
#include <chrono>

// generate a single instruction, unlike __builtin_readcyclecounter()
inline long rdtsc() {
    long ret;
    asm volatile ("rdtsc" : "=A"(ret));
    return ret;
}

struct Timestamp {
	Timestamp() : m_ticks(rdtsc()) { }

	long elapsed(Timestamp a = Timestamp()) const { return a.m_ticks - m_ticks; }

	double elapsed_ms(Timestamp a = Timestamp()) const { return elapsed(a) * s_milisec_per_tick; }

	static void init();

private:
	long m_ticks;
	static double s_milisec_per_tick;
};
