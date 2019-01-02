#include <catch.hpp>
#include <geom/broadphase.h>
#include <core/range.h>
#include <core/timestamp.h>

void Bruteforce(const vector<sphere>& spheres, vector<pair<int, int>>& intersections) {
	for (auto i : range(spheres.size()))
		for (auto j : range(i)) {
			const sphere& a = spheres[i];
			const sphere& b = spheres[j];
			if (squared(a.center() - b.center()) <= squared(a.radius() + b.radius()))
				intersections.emplace_back(i, j);
		}
}

TEST_CASE("HashBroadphase", "[.][broadphase][benchmark]") {
	HashBroadphase bp;
	std::default_random_engine rnd;
	for (int t : range(10)) {
		vector<sphere> spheres;
		for (int i : range(20000))
			spheres.push_back(sphere(uniform3p(rnd, -10, 10), 0.02));
		bp.reset(1.05 * 0.02 / 0.3 * 8, 0.02);

		vector<pair<int, int>> result1;
		vector<pair<int, int>> result2;
		vector<int> result;

		Timestamp ta;
		Bruteforce(spheres, result1);
		Timestamp tb;
		for (int i : range(spheres.size())) {
			const sphere& a = spheres[i];
			bp.getCandidates(a, result);
			for (int j : result) {
				const sphere& b = spheres[j];
				if (squared(a.center() - b.center()) <= squared(a.radius() + b.radius()))
					result2.emplace_back(i, j);
			}
			result.clear();
			bp.add(a, i);
		}
		Timestamp tc;

		print("Brute %.1f %s Hash %.1f %s\n", ta.elapsed_ms(tb), result1.size(), tb.elapsed_ms(tc), result2.size());
	}
}
