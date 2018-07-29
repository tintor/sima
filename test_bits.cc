#include "bits.h"
#include "catch.hpp"
#include <iostream>

using namespace std::literals;
using namespace std;

TEST_CASE("is_power2") {
	for (uint i = 0; i <= 2000; i++) {
		if (i <= 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
			REQUIRE(is_power2(i));
		else
			REQUIRE(!is_power2(i));
	}
}

TEST_CASE("round_up_power2") {
	REQUIRE(0 == round_up_power2(0lu));
	REQUIRE(1 == round_up_power2(1lu));
	REQUIRE(2 == round_up_power2(2lu));
	REQUIRE(4 == round_up_power2(3lu));
	REQUIRE(4 == round_up_power2(4lu));
	REQUIRE(8 == round_up_power2(5lu));

	for (ulong i = 1; i < 100000; i++) {
		ulong a = round_up_power2(i);
		ulong b = round_up_power2(i + 1);
		if (a != b)
			REQUIRE(a + a == b);
	}

	for (ulong i = 2; i <= 63; i++) {
		ulong a = 1;
		a <<= i;
		REQUIRE(round_up_power2(a - 1) == a);
		REQUIRE(round_up_power2(a) == a);
		REQUIRE(round_up_power2(a + 1) == a + a);
	}
}

TEST_CASE("clz(uint)") {
	uint a = 0x1;
	for (int i = 0; i <= 31; i++)
		REQUIRE(clz(a << i) == 31 - i);
}

TEST_CASE("clz(ulong)") {
	ulong a = 0x1;
	for (int i = 0; i <= 63; i++)
		REQUIRE(clz(a << i) == 63 - i);
}

TEST_CASE("clz(uint128)") {
	uint128 a = 0x1;
	for (int i = 0; i <= 127; i++)
		REQUIRE(clz(a << i) == 127 - i);
}

TEST_CASE("popcount(uint)") {
	uint a = 0x1;
	REQUIRE(popcount(a ^ a) == 0);
	REQUIRE(popcount(a) == 1);
	REQUIRE(popcount((a << 31) | 1) == 2);
}

TEST_CASE("popcount(ulong)") {
	ulong a = 0x1;
	REQUIRE(popcount(a ^ a) == 0);
	REQUIRE(popcount(a) == 1);
	REQUIRE(popcount((a << 63) | 1) == 2);
}

TEST_CASE("popcount(uint128)") {
	uint128 a = 0x1;
	REQUIRE(popcount(a ^ a) == 0);
	REQUIRE(popcount(a) == 1);
	REQUIRE(popcount((a << 127) | 1) == 2);
}

TEST_CASE("Bits() capacity() data() size() set() reset()") {
	Bits b;
	REQUIRE(b.capacity() == 0);
	REQUIRE(b.data() == nullptr);

	b.set(0);
	b.set(2);
	REQUIRE(b.capacity() == 1);
	REQUIRE(b.data()[0] == 5u);
	REQUIRE(clz(5lu) == 61);
	REQUIRE(b.size() == 3);

	b.set(64);
	REQUIRE(b.capacity() == 2);
	REQUIRE(b.size() == 65);
	b.reset(640);
	REQUIRE(b.capacity() == 2);
	b.reset(64);
	REQUIRE(b.capacity() == 2);
	REQUIRE(b.data()[0] == 5u);
	REQUIRE(b.data()[1] == 0u);
	REQUIRE(b.size() == 3);
}

TEST_CASE("Bits(const Bits&) size_words() get_word() set_word()") {
	Bits b;
	REQUIRE(b.size_words() == 0);
	b.set_word(0, 1);
	b.set_word(0, 0);
	REQUIRE(b.capacity() == 1);
	REQUIRE(b.size_words() == 0);

	b.set_word(0, 23);
	REQUIRE(b.get_word(0) == 23);
	REQUIRE(b.get_word(1) == 0);
	b.set_word(1, 17);
	REQUIRE(b.capacity() == 2);
	REQUIRE(b.data()[0] == 23);
	REQUIRE(b.data()[1] == 17);
	REQUIRE(b.size_words() == 2);

	b.set_word(1, 0);
	REQUIRE(b.capacity() == 2);
	REQUIRE(b.data()[0] == 23);
	REQUIRE(b.data()[1] == 0);
	REQUIRE(b.size_words() == 1);

	b.set_word(2, 0);
	REQUIRE(b.capacity() == 2);
	REQUIRE(b.data()[0] == 23);
	REQUIRE(b.data()[1] == 0);

	Bits c(b);
	REQUIRE(c.capacity() == 1);
	REQUIRE(c.data()[0] == 23);
}

TEST_CASE("popcount()") {
	Bits b;
	REQUIRE(b.popcount() == 0);
	b.set_word(1, 5);
	b.set_word(2, 5);
	b.set_word(2, 0);
	REQUIRE(b.popcount() == 2);
}

TEST_CASE("Bits(uint32_t) empty") {
	uint32_t a = 0;
	Bits p(a);
	REQUIRE(p.capacity() == 0);
	REQUIRE(p.data() == nullptr);
}

TEST_CASE("Bits(uint32_t)") {
	uint32_t a = 0x12345678;
	Bits p(a);
	REQUIRE(p.capacity() == 1);
	REQUIRE(p.data() != nullptr);
	REQUIRE(p.data()[0] == a);
}

TEST_CASE("Bits(uint64_t) empty") {
	ulong a = 0;
	Bits p(a);
	REQUIRE(p.capacity() == 0);
	REQUIRE(p.data() == nullptr);
}

TEST_CASE("Bits(uint64_t)") {
	ulong a = 0x1234567878654321lu;
	Bits p(a);
	REQUIRE(p.capacity() == 1);
	REQUIRE(p.data() != nullptr);
	REQUIRE(p.data()[0] == a);
}

TEST_CASE("Bits(uint128_t) empty") {
	ulong a = 0;
	Bits p(a);
	REQUIRE(p.capacity() == 0);
	REQUIRE(p.data() == nullptr);
}

TEST_CASE("Bits(uint128_t) small") {
	ulong a = 0x1234567878654321lu;
	uint128 A = a;
	Bits p(A);
	REQUIRE(p.capacity() == 1);
	REQUIRE(p.data() != nullptr);
	REQUIRE(p.data()[0] == a);
}

TEST_CASE("Bits(uint128_t) big") {
	uint64_t a = 3;
	uint64_t b = 7;
	__uint128_t A = (static_cast<__uint128_t>(a) << 64) | b;
	REQUIRE(clz(A) == 62);
	Bits p(A);
	REQUIRE(p.capacity() == 2);
	REQUIRE(p.data() != nullptr);
	REQUIRE(p.data()[1] == a);
	REQUIRE(p.data()[0] == b);
}

TEST_CASE("Bits(string_view)") {
	REQUIRE(Bits(""sv).capacity() == 0);
	REQUIRE(Bits("0"sv).capacity() == 0);
	REQUIRE(Bits("00"sv).capacity() == 0);
	REQUIRE(Bits("000"sv).capacity() == 0);

	REQUIRE(Bits("7"sv).capacity() == 1);
	REQUIRE(Bits("7"sv).get_word(0) == 7);

	REQUIRE(Bits("4a"sv).capacity() == 1);
	REQUIRE(Bits("4a"sv).get_word(0) == 0x4a);

	REQUIRE(Bits("4a2"sv).capacity() == 1);
	REQUIRE(Bits("4a2"sv).get_word(0) == 0x4a2);

	std::string_view s = "f000000000000000a"sv;
	REQUIRE(Bits(s).capacity() == 2);
	REQUIRE(Bits(s).get_word(1) == 0xf);
	REQUIRE(Bits(s).get_word(0) == 0xa);
}
