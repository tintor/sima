#pragma once
#include <string>
#include <memory>
#include <array>

class Callstack {
public:
	Callstack();
	static std::unique_ptr<char, void(*)(char*)> demangle(const char* symbol);
	void write(std::string& out) const;

private:
	std::array<void*, 20> _stack;
	int _size;
};
