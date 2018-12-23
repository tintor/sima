#include <lisp.h>
#include <core/format.h>
#include <core/string_util.h>
#include <catch.hpp>

void test(bool identical, string_view cases) {
	for (auto line : split(cases, '\n')) {
		if (line.size() == 0)
			continue;
		auto parts = split(line, '\t');
		Lisp lisp;

		auto ap = lisp.parse(parts[0]);
		auto a = lisp.eval(ap);
		auto bp = lisp.parse(parts.back());
		auto b = lisp.eval(bp);

		bool result = identical ? a == b : Lisp::equal(a, b);
		if (!result) {
			print("as = %s\n", parts[0]);
			print("bs = %s\n", parts.back());
			print("ap = %s\n", lisp.format(ap));
			print("bp = %s\n", lisp.format(bp));
			print("a = %s\n", lisp.format(a));
			print("b = %s\n", lisp.format(b));
			print("a not %s to b\n", identical ? "identical" : "equal");
		}
		REQUIRE(result);
	}
}

TEST_CASE("lisp identical", "[lisp]") {
test(true, R"""(
(# (+ 2 1) This is comment)		()

0			0
1			1
-1			-1
true		true
false		false
()			()
(+ 1 2)		3
(- 5 1)		4
(+ -7 -10)	-17
(- 2)		-2
(- -2)		2

(' a)		(' a)
(' if)		(' if)

(> 1 0)		true
(> 0 1)		false

(== 0 0)		true
(== 1 1)		true
(== 0 1)		false
(== 1 0)		false
(== 0 "e")		false
(== "e" 0)		false
(== "e" "e")	true
(== "e" "a")	false

(== () ())		true
(== () (1))		false
(== (1) (1))	true
(== (1 2) (1 3))	false
(== (1 2) (1 2))	true

(if true 2 3)		2
(if false 2 3)		3
(if (> 1 0) 2 3)	2
(if (> 0 1) 2 3)	3

(if false 3)		()

(not true)		false
(not false)		true
(not (> 1 0))	false
(not (> 0 1))	true

(and true true)		true
(and true false)	false
(and false false)	false
(and false false)	false

(and true true true true)	true
(and false true true true)	false
(and true true true false)	false

(and true (> 1 0))	true
(and true (> 0 1))	false

((= a (+ 2 3)) a)	5

(parse "a")				(' a)

(eval (parse "(+ 1 2)"))		3

((= a 3) (eval (parse "a")))	3

((= i 2) (+= i 1))		3
((= i 0) (while false (= i 1)) i)		0
((= i 0) (while (< i 10) (+= i 1)) i)	10

(while true break)		()

(=== 2 2)			true
(=== "ma" "ma")		false
(=== (' a) (' a))	true
(=== () ())			true
(=== (1 2) (1 2))	false
(=== true true)		true
(=== false false)	true
(=== true false)	false
(=== false true)	false

(size ())		0
(size (3))		1
(size (5 6))	2
(size "")		0
(size "a")		1
(size "ab")		2
(size (' z))	1
(size (' xyz))	3
(size (parse "(4 5 6 7)"))	4

(size (map))				0
(size (map "a" 1))			1
(size (map "a" 1 "b" 2))	2

(in 10 (map 10 1))		true
(in 20 (map 10 1))		false

(in "a" (map "a" 1))	true
(in "b" (map "a" 1))	false

(get (10 20 30) -1)		()
(get (10 20 30) 0)		10
(get (10 20 30) 1)		20
(get (10 20 30) 2)		30
(get (10 20 30) 3)		()

(get (map "a" 1) "a")	1
(get (map "a" 1) "b")	()

((= m (map)) (set m "a" 4) (set m "b" 3) (get m "a"))			4
((= m (map)) (set m "a" 4) (set m "a" 3) (+ (get m "a") 10))	13
)""");
}

TEST_CASE("lisp equal", "[lisp]") {
	test(false, R"""(
(> "e" 0)		(error "(>) arg not int" "e")
(not 0)			(error "(not) arg not bool" 0)
(and 0 true)	(error "(and) arg not bool" 0)
(or 0 true)		(error "(or) arg not bool" 0)

(parse "(+ 1 2)")		(' (+ 1 2))
)""");
}

static string eval(string_view s) {
	Lisp lisp;
	auto a = lisp.eval(lisp.parse(s));
	return lisp.format(a);
}

TEST_CASE("lisp format", "[lisp]") {
	REQUIRE(eval("(map)") == "(map)");
	REQUIRE(eval("(map 7 10)") == "(map 7 10)");
	auto s = eval("(map 3 4 () -2)");
	REQUIRE((s == "(map 3 4 () -2)" || s == "(map () -2 3 4)"));

	REQUIRE(eval("((= m (1 2 3)) (set m 0 10) m)") == "(10 2 3)");
	REQUIRE(eval("((= m (1 2 3)) (set m -1 40) m)") == "(error \"(set) index out of bounds\" -1)");
}
