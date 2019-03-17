#include <catch.hpp>
#include <core/format.h>
#include <core/bplustree.h>

TEST_CASE("bplustree3 random inserts", "[bplustree]") {
	std::default_random_engine rnd;
	for (auto x : range(2500)) {
		bplustree_map<int, int, 3> b;
		REQUIRE(b.size() == 0);
		map<int, int> e;
		std::uniform_int_distribution dist(0, x-1);
		for (auto i : range(x)) {
			int k = dist(rnd);
			b[k] = -k;
			e[k] = -k;
			REQUIRE(b.size() == e.size());
		}

		auto bi = b.begin();
		auto ei = e.begin();
		uint i = 0;
		while (ei != e.end()) {
			REQUIRE(bi != b.end());
			ASSERT_ALWAYS(*bi == *ei, "i=%s b=%s e=%s", i, bi->first, ei->first);
			++bi;
			++ei;
			i += 1;
		}
		REQUIRE(bi == b.end());
	}
}

TEST_CASE("bplustree4 random inserts", "[bplustree]") {
	std::default_random_engine rnd;
	for (auto x : range(2500)) {
		bplustree_map<int, int, 4> b;
		REQUIRE(b.size() == 0);
		map<int, int> e;
		std::uniform_int_distribution dist(0, x-1);
		for (auto i : range(x)) {
			int k = dist(rnd);
			b[k] = -k;
			e[k] = -k;
			REQUIRE(b.size() == e.size());
		}

		auto bi = b.begin();
		auto ei = e.begin();
		uint i = 0;
		while (ei != e.end()) {
			REQUIRE(bi != b.end());
			ASSERT_ALWAYS(*bi == *ei, "i=%s b=%s e=%s", i, bi->first, ei->first);
			++bi;
			++ei;
			i += 1;
		}
		REQUIRE(bi == b.end());
	}
}

constexpr uint C = 10 * 1000 * 1000;

template<typename T>
void test() {
	std::default_random_engine rnd;
	rnd.seed(0);
	T b;
	std::uniform_int_distribution dist(0u, C-1);
	for (auto i : range(C)) {
		int k = dist(rnd);
		b[k] = -k;
	}

	/*uint sum = 0;
	auto it = b.begin();
	while (it != b.end()) {
		sum += it->second;
		++it;
	}
	print("%s\n", sum);*/
}

TEST_CASE("map random inserts - benchmark", "[.][bplustree]") {
	test<map<int, int>>();
}

TEST_CASE("unordered_map random inserts - benchmark", "[.][bplustree]") {
	test<unordered_map<int, int>>();
}

TEST_CASE("b+tree8 random inserts - benchmark", "[.][bplustree]") {
	test<bplustree_map<int, int, 8>>();
}

TEST_CASE("b+tree12 random inserts - benchmark", "[.][bplustree]") {
	test<bplustree_map<int, int, 12>>();
}

TEST_CASE("b+tree16 random inserts - benchmark", "[.][bplustree]") {
	test<bplustree_map<int, int, 16>>();
}

TEST_CASE("b+tree20 random inserts - benchmark", "[.][bplustree]") {
	test<bplustree_map<int, int, 20>>();
}

TEST_CASE("b+tree24 random inserts - benchmark", "[.][bplustree]") {
	test<bplustree_map<int, int, 24>>();
}

TEST_CASE("b+tree28 random inserts - benchmark", "[.][bplustree]") {
	test<bplustree_map<int, int, 28>>();
}

TEST_CASE("b+tree32 random inserts - benchmark", "[.][bplustree]") {
	test<bplustree_map<int, int, 32>>();
}

TEST_CASE("ranked_set", "[ranked_set]") {
	std::default_random_engine rnd;
	std::array<uint, 1000> a;
	for (auto i : range(a.size()))
		a[i] = i * 10;
	std::shuffle(a.begin(), a.end(), rnd);
	ranked_set<int, 5> b;
	for (auto e : a)
		b.insert(e);
	for (auto i : range(a.size()))
		print("%s %s\n", i, b.findByRank(i));
}
