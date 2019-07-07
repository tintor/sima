#pragma once

struct Timestamp {
	Timestamp() : _ticks(__builtin_readcyclecounter()) { }
	Timestamp(long a) : _ticks(a) { }

	long elapsed(Timestamp a = Timestamp()) const { return a._ticks - _ticks; }

	double elapsed_s(Timestamp a = Timestamp()) const { return elapsed(a) * _ms_per_tick * 1e-3; }
	double elapsed_ms(Timestamp a = Timestamp()) const { return elapsed(a) * _ms_per_tick; }

	static void init();

	static double to_s(long ticks) { return ticks * _ms_per_tick * 1e-3; }

private:
	long _ticks;
	static double _ms_per_tick;
};
