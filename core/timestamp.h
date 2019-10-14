#pragma once
#include <core/std.h>
#include <core/format.h>

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

class AutoTimestamp {
public:
	AutoTimestamp(string_view name) : _name(name) {}
	~AutoTimestamp() { print("%s took %s ms\n", _name, _begin.elapsed_ms()); }
private:
	string_view _name;
	Timestamp _begin;
};

template<typename Func>
auto measure(string_view name, Func func) {
	AutoTimestamp mm(name);
	return func();
}

#define MEASURE(X) measure(#X, [&](){ return X; } )
