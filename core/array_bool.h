#pragma once
#include <array>

template<uint Size>
struct array_bool {
	static constexpr uint Words = (Size + 31) / 32;
	std::array<uint, Words> words;

	static uint size() {
		return Size;
	}

	bool operator[](int index) const {
		uint mask = uint(1) << (index % 32);
		return (words[index / 32] & mask) != 0;
	}

	struct BitRef {
		uint mask;
		uint* word;

		operator bool() {
			return (*word & mask) != 0;
		}

		BitRef& operator=(bool a) {
			if (a)
				*word |= mask;
			else
				*word &= ~mask;
			return *this;
		}
	};

	BitRef operator[](int index) {
		return { uint(1) << (index % 32), words.data() + index / 32 };
	}

	void set(int index) {
		uint mask = uint(1) << (index % 32);
		words[index / 32] |= mask;
	}

	void reset(int index) {
		uint mask = uint(1) << (index % 32);
		words[index / 32] &= ~mask;
	}
};

