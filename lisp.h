#pragma once
#include <core/std.h>

class Lisp {
public:
	using Node = const void*;

	Node parse(string_view code);
	Node eval(Node n);

	Node parse_block(string_view code);
	Node eval_block(Node n);

	void format(string& s, Node node);

	bool equal(Node, Node);

private:
	void* m_handle;
};
