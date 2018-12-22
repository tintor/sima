#include <lisp.h>
#include <core/range.h>
#include <core/format.h>
#include <core/exception.h>

static bool is_cent(string_view s) {
	if (s.size() > 0 && s[0] == '-')
		s = s.substr(1);
	if (s.empty())
		return false;
	for (char c : s)
		if (c < '0' || c > '9')
			return false;
	return true;
}

static cent parse_cent(string_view s) {
	cent a = 0;
	int i = 0;
	if (s[0] == '-') {
		i += 1;
		while (i < s.size()) {
			char c = s[i];
			a = a * 10 + c - '0';
			i += 1;
		}
		a = -a;
	} else {
		while (i < s.size()) {
			char c = s[i];
			a = a * 10 + c - '0';
			i += 1;
		}
	}
	return a;
}

struct Node;
static unordered_map<string_view, const Node*> g_symbols;

enum class Type : uint8_t {
	None=0, List, Keyword, Symbol, Integer, Decimal, String
};

struct Node {
private:
	Type _type = Type::None;
	uint _size = 0;
	void* _data = nullptr; // TODO allocate memory inline

	static Node* alloc_raw(uint size) {
		Node* n = new Node;
		n->_size = size;
		n->_data = malloc(size);
		return n;
	}

	static Node* alloc_list(uint size) {
		Node* n = new Node;
		n->_type = Type::List;
		n->_size = size;
		n->_data = malloc(size * 8);
		return n;
	}

	const Node* get(uint index) const {
		ASSERT_ALWAYS(is_list());
		ASSERT_ALWAYS(index < _size);
		const Node** data = (const Node**) _data;
		return data[index];
	}

	void set(uint index, const Node* n) {
		ASSERT_ALWAYS(is_list());
		ASSERT_ALWAYS(index < _size);
		const Node** data = (const Node**) _data;
		data[index] = n;
	}

public:
	static const Node* symbol(string_view a, bool reserved = true) {
		auto it = g_symbols.find(a);
		if (it != g_symbols.end())
			return it->second;

		Node* n = alloc_raw(a.size());
		n->_type = reserved ? Type::Keyword : Type::Symbol;
		memcpy(n->_data, a.data(), a.size());
		g_symbols.insert({n->symbol(), n});
		return n;
	}

	static const Node* integer(cent a) {
		Node* n = alloc_raw(sizeof(a));
		n->_type = Type::Integer;
		memcpy(n->_data, &a, sizeof(a));
		return n;
	}

	static const Node* str(string_view a) {
		Node* n = alloc_raw(a.size());
		n->_type = Type::String;
		memcpy(n->_data, a.data(), a.size());
		return n;
	}

	static const Node* list() { return alloc_list(0); }

	static const Node* list(string_view a, string_view b) {
		auto n = alloc_list(2);
		n->set(0, symbol(a));
		n->set(1, str(b));
		return n;
	}

	static const Node* list(string_view a, string_view b, const Node* c) {
		auto n = alloc_list(3);
		n->set(0, symbol(a));
		n->set(1, str(b));
		n->set(2, c);
		return n;
	}

	static const Node* list(const Node* a, string_view b) {
		auto n = alloc_list(2);
		n->set(0, a);
		n->set(1, str(b));
		return n;
	}

	static const Node* list(const Node* a, const Node* b) {
		auto n = alloc_list(2);
		n->set(0, a);
		n->set(1, b);
		return n;
	}

	static const Node* list(const Node* a, string_view b, const Node* c) {
		auto n = alloc_list(3);
		n->set(0, a);
		n->set(1, str(b));
		n->set(2, c);
		return n;
	}

	static const Node* list(span<const Node*> elements) {
		auto n = alloc_list(elements.size());
		for (uint i : range(elements.size()))
			n->set(i, elements[i]);
		return n;
	}

	static const Node* list(std::initializer_list<const Node*> elements) {
		auto n = alloc_list(elements.size());
		for (uint i : range(elements.size()))
			n->set(i, elements.begin()[i]);
		return n;
	}

	bool is_symbol() const {
		return _type == Type::Keyword || _type == Type::Symbol;
	}
	bool is_strict_symbol() const {
		return _type == Type::Symbol;
	}
	bool is_keyword() const {
		return _type == Type::Keyword;
	}
	string_view symbol() const {
		ASSERT_ALWAYS(is_strict_symbol() || is_keyword());
		const char* atom = (const char*) _data;
		return string_view(atom, _size);
	}

	bool is_integer() const {
		return _type == Type::Integer;
	}
	cent integer() const {
		ASSERT_ALWAYS(is_integer());
		return *(cent*)_data;
	}

	bool is_string() const {
		return _type == Type::String;
	}
	string_view str() const {
		ASSERT_ALWAYS(is_string());
		const char* atom = (const char*) _data;
		return string_view(atom, _size);
	}

	bool is_list() const {
		return _type == Type::List;
	}
	uint size() const {
		ASSERT_ALWAYS(is_list());
		return _size;
	}
	const Node* first() const {
		ASSERT_ALWAYS(is_list());
		ASSERT_ALWAYS(size() > 0);
		const Node** data = (const Node**) _data;
		return data[0];
	}
	auto elements() const {
		ASSERT_ALWAYS(is_list());
		const Node** data = (const Node**) _data;
		return span<const Node*>(data, _size);
	}
	auto args() const {
		ASSERT_ALWAYS(size() > 0);
		return elements().pop_front();
	}
};

const Node* NUMBER = Node::symbol("#");

const Node* ERROR = Node::symbol("error");
const Node* DEF = Node::symbol("def");

const Node* TRUE = Node::symbol("true");
const Node* FALSE = Node::symbol("false");
const Node* AND = Node::symbol("and");
const Node* OR = Node::symbol("or");
const Node* NOT = Node::symbol("not");

const Node* AMP = Node::symbol("&");
const Node* PIPE = Node::symbol("|");
const Node* CAPPA = Node::symbol("^");
const Node* TILDA = Node::symbol("~");

const Node* PLUS = Node::symbol("+");
const Node* PLUS_EQUAL = Node::symbol("+=");
const Node* MINUS = Node::symbol("-");
const Node* MINUS_EQUAL = Node::symbol("-=");
const Node* STAR = Node::symbol("*");
const Node* STAR_EQUAL = Node::symbol("*=");
const Node* DIV_DIV = Node::symbol("//");
const Node* DIV_DIV_EQUAL = Node::symbol("//=");
const Node* DIV = Node::symbol("/");
const Node* DIV_EQUAL = Node::symbol("/=");
const Node* PERCENT = Node::symbol("%");
const Node* PERCENT_EQUAL = Node::symbol("%=");

const Node* LESS_LESS = Node::symbol("<<");
const Node* GREATER_GREATER = Node::symbol(">>");

const Node* IF = Node::symbol("if");
const Node* WHILE = Node::symbol("while");
const Node* BREAK = Node::symbol("break");
const Node* CONTINUE = Node::symbol("continue");

const Node* LESS = Node::symbol("<");
const Node* LESS_EQUAL = Node::symbol("<=");
const Node* GREATER = Node::symbol(">");
const Node* GREATER_EQUAL = Node::symbol(">=");
const Node* EQUAL_EQUAL = Node::symbol("==");
const Node* NOT_EQUAL = Node::symbol("!=");

const Node* PRINT = Node::symbol("print");
const Node* PRINTLN = Node::symbol("println");
const Node* EQUAL = Node::symbol("=");

const Node* N_EMPTY = Node::list(); // TODO make all empty lists return the same pointer
const Node* N_BREAK = Node::list(ERROR, BREAK);
const Node* N_CONTINUE = Node::list(ERROR, CONTINUE);

static bool space(char c) { return c == ' ' || c == '\t' || c == '\n'; }

static const Node* parse(string_view code) {
	int p = 0;
	vector<const Node*> stack = {};
	vector<uint> size = {0};
	while (p < code.size()) {
		char c = code[p];
		if (space(c)) {
			p += 1;
			continue;
		}
		if (c == '(') {
			size.push_back(0);
			p += 1;
			continue;
		}
		if (c == ')') {
			if (size.size() == 1) {
				print("unexpected ')'\n");
				return nullptr;
			}
			uint s = size.back();
			size.pop_back();
			const Node* n = Node::list(span<const Node*>(stack.data() + stack.size() - s, s));
			stack.resize(stack.size() - s);

			stack.push_back(n);
			size.back() += 1;
			p += 1;
			continue;
		}
		if (c == '"') {
			p += 1;
			int s = p;
			while (p < code.size() && code[p] != '"')
				p += 1;
			if (p == code.size()) {
				print("unterminated quote\n");
				return nullptr;
			}
			const Node* n = Node::str(code.substr(s, p - s));
			stack.push_back(n);
			size.back() += 1;
			p += 1;
			continue;
		}
		// anything else is atom
		int s = p;
		while (p < code.size() && !space(code[p]) && code[p] != '(' && code[p] != ')')
			p += 1;
		string_view ss = code.substr(s, p - s);
		const Node* n = is_cent(ss) ? Node::integer(parse_cent(ss)) : Node::symbol(ss, false);
		stack.push_back(n);
		size.back() += 1;
	}
	if (size.size() != 1) {
		print("unterminated list\n");
		return nullptr;
	}
	return Node::list(span<const Node*>(stack.data(), stack.size()));
}

int length(const Node* n) {
	ASSERT_ALWAYS(n != nullptr);
	if (n->is_integer()) {
		int s = 0;
		cent c = n->integer();
		if (c <= 0) {
			s += 1;
			c = -c;
		}
		while (c > 0) {
			s += 1;
			c /= 10;
		}
		return s;
	}
	if (n->is_string()) {
		int s = 0;
		for (char c : n->str())
			if (c == '"' || c == '\\')
				s += 1;
		return s + 2 + n->str().size();
	}
	if (n->is_strict_symbol() || n->is_keyword())
		return n->symbol().size();

	ASSERT_ALWAYS(n->is_list());
	if (n->size() == 0)
		return 2;
	int s = 0;
	for (const Node* e : n->elements())
		s += length(e);
	return s + 1 + n->size();
}

void format_e(string& s, string_view spec, const Node* n) {
	ASSERT_ALWAYS(n != nullptr);
	if (n->is_integer()) {
		format_e(s, "%s", n->integer());
		return;
	}
	if (n->is_strict_symbol() || n->is_keyword()) {
		s += n->symbol();
		return;
	}
	if (n->is_string()) {
		s += '"';
		for (char c : n->str()) {
			if (c == '"' || c == '\\')
				s += '\\';
			s += c;
		}
		s += '"';
		return;
	}
	ASSERT_ALWAYS(n->is_list());

	if (length(n) <= 120) {
		s += '(';
		bool first = true;
		for (const Node* e : n->elements()) {
			if (first)
				first = false;
			else
				s += ' ';
			format_e(s, spec, e);
		}
		s += ')';
		return;
	}

	static int s_depth = 0;
	for (int i : range(s_depth * 2))
		s += ' ';
	s += "(\n";
	s_depth += 1;

	for (const Node* e : n->elements()) {
		ASSERT_ALWAYS(e != nullptr);
		if (!e->is_list() || length(e) <= 120) {
			for (int i : range(s_depth * 2))
				s += ' ';
			format_e(s, spec, e);
		} else {
			format_e(s, spec, e);
		}
		s += '\n';
	}

	s_depth -= 1;
	for (int i : range(s_depth * 2))
		s += ' ';
	s += ")";
}

static unordered_map<const Node*, const Node*> g_variables;

static bool is_error(const Node* e) {
	return e->is_list() && e->size() > 0 && e->elements()[0] == ERROR;
}

static const Node* eval(const Node* n);

#define FAIL(...) return Node::list(ERROR, ##__VA_ARGS__)

const Node* eval_plus(span<const Node*> args) {
	if (args.size() < 2)
		FAIL("(+) min 2 args");
	cent acc = 0;
	for (const Node* e : args) {
		e = eval(e);
		if (is_error(e))
			return e;
		if (!e->is_integer())
			FAIL("(+) arg not int", e);
		acc += e->integer();
	}
	return Node::integer(acc);
}

const Node* eval_plus_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(+=) exactly 2 args");

	const Node* a = args[0];
	if (!a->is_strict_symbol())
		FAIL("(+=) lhs must be symbol", a);
	auto it = g_variables.find(a);
	if (it == g_variables.end())
		FAIL("(+=) symbol not defined", a);
	auto& av = it->second;
	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(+=) arg not int", b);
	b = Node::integer(av->integer() + b->integer());
	av = b;
	return b;
}

const Node* eval_minus(span<const Node*> args) {
	if (args.size() == 1) {
		auto a = eval(args[0]);
		if (is_error(a))
			return a;
		if (!a->is_integer())
			FAIL("(-) arg not int", a);
		return Node::integer(-a->integer());
	}

	if (args.size() == 2) {
		auto a = eval(args[0]);
		if (is_error(a))
			return a;
		if (!a->is_integer())
			FAIL("(-) arg not int", a);

		auto b = eval(args[1]);
		if (is_error(b))
			return b;
		if (!b->is_integer())
			FAIL("(-) arg not int", b);

		return Node::integer(a->integer() - b->integer());
	}
	FAIL("(-) exactly 1 or 2 args");
}

const Node* eval_minus_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(-=) exactly 2 args");

	const Node* a = args[0];
	if (!a->is_strict_symbol())
		FAIL("(-=) lhs must be symbol", a);
	auto it = g_variables.find(a);
	if (it == g_variables.end())
		FAIL("(+=) symbol not defined", a);
	auto& av = it->second;
	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(-=) arg not int", b);
	b = Node::integer(av->integer() - b->integer());
	av = b;
	return b;
}

const Node* eval_not(span<const Node*> args) {
	if (args.size() != 1)
		FAIL("(not) exactly 1 arg");

	auto a = eval(args[0]);
	if (is_error(a))
		return a;
	if (a == TRUE)
		return FALSE;
	if (a == FALSE)
		return TRUE;
	FAIL("(not) arg not bool", a);
}

const Node* eval_and(span<const Node*> args) {
	if (args.size() == 0)
		FAIL("(and) min 1 arg");

	for (auto e : args) {
		e = eval(e);
		if (is_error(e))
			return e;
		if (e == FALSE)
			return FALSE;
		if (e != TRUE)
			FAIL("(and) arg not bool", e);
	}
	return TRUE;
}

const Node* eval_or(span<const Node*> args) {
	if (args.size() == 0)
		FAIL("(or) min 1 arg");

	for (auto e : args) {
		e = eval(e);
		if (is_error(e))
			return e;
		if (e == TRUE)
			return TRUE;
		if (e != FALSE)
			FAIL("(or) arg not bool", e);
	}
	return FALSE;
}

const Node* eval_star(span<const Node*> args) {
	if (args.size() < 2)
		FAIL("(*) min 2 args");
	cent acc = 1;
	for (auto e : args) {
		e = eval(e);
		if (is_error(e))
			return e;
		if (!e->is_integer())
			FAIL("(*) arg not int", e);
		acc *= e->integer();
	}
	return Node::integer(acc);
}

const Node* eval_star_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(*=) exactly 2 args");

	const Node* a = args[0];
	if (!a->is_strict_symbol())
		FAIL("(*=) lhs must be symbol", a);
	auto it = g_variables.find(a);
	if (it == g_variables.end())
		FAIL("(*=) symbol not defined", a);
	auto& av = it->second;
	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(*=) arg not int", b);
	b = Node::integer(av->integer() * b->integer());
	av = b;
	return b;
}

const Node* eval_div_div(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(//) exactly 2 args");

	auto a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(//) arg not int", a);

	auto b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(//) arg not int", b);

	return Node::integer(a->integer() / b->integer());
}

const Node* eval_div_div_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(//=) exactly 2 args");

	const Node* a = args[0];
	if (!a->is_strict_symbol())
		FAIL("(//=) lhs must be symbol", a);
	auto it = g_variables.find(a);
	if (it == g_variables.end())
		FAIL("(//=) symbol not defined", a);
	auto& av = it->second;
	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(//=) arg not int", b);
	b = Node::integer(av->integer() / b->integer());
	av = b;
	return b;
}

const Node* eval_percent(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(%) exactly 2 args");

	auto a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(%) arg not int", a);

	auto b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(%) arg not int", b);

	return Node::integer(a->integer() % b->integer());
}

const Node* eval_percent_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(%=) exactly 2 args");

	const Node* a = args[0];
	if (!a->is_strict_symbol())
		FAIL("(%=) lhs must be symbol", a);
	auto it = g_variables.find(a);
	if (it == g_variables.end())
		FAIL("(%=) symbol not defined", a);
	auto& av = it->second;
	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(%=) arg not int", b);
	b = Node::integer(av->integer() % b->integer());
	av = b;
	return b;
}

const Node* eval_less_less(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(<<) exactly 2 args");

	auto a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(<<) arg not int", a);

	auto b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(<<) arg not int", b);

	return Node::integer(a->integer() << b->integer());
}

const Node* eval_greater_greater(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(>>) exactly 2 args");

	auto a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(>>) arg not int", a);

	auto b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(>>) arg not int", b);

	return Node::integer(a->integer() >> b->integer());
}

const Node* eval_less(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(<) exactly 2 args");

	const Node* a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(<) arg not int", a);

	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(<) arg not int", b);

	return (a->integer() < b->integer()) ? TRUE : FALSE;
}

const Node* eval_less_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(<=) exactly 2 args");

	const Node* a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(<=) arg not int", a);

	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(<=) arg not int", b);

	return (a->integer() <= b->integer()) ? TRUE : FALSE;
}

const Node* eval_greater(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(>) exactly 2 args");

	const Node* a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(>) arg not int", a);

	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(>) arg not int", b);

	return (a->integer() > b->integer()) ? TRUE : FALSE;
}

const Node* eval_greater_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(>=) exactly 2 args");

	const Node* a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(>=) arg not int", a);

	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(>=) arg not int", b);

	return (a->integer() >= b->integer()) ? TRUE : FALSE;
}

// TODO compare any type (string, decimal, bool)!
const Node* eval_equal_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(==) exactly 2 args");

	const Node* a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(==) arg not int", a);

	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(==) arg not int", b);

	return (a->integer() == b->integer()) ? TRUE : FALSE;
}

const Node* eval_not_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(!=) exactly 2 args");

	const Node* a = eval(args[0]);
	if (is_error(a))
		return a;
	if (!a->is_integer())
		FAIL("(!=) arg not int", a);

	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	if (!b->is_integer())
		FAIL("(!=) arg not int", b);

	return (a->integer() != b->integer()) ? TRUE : FALSE;
}

const Node* eval_print(span<const Node*> args, bool endl) {
	for (auto e : args) {
		e = eval(e);
		if (is_error(e))
			return e;
		if (e->is_string())
			print("%s", e->str());
		else
			print("%s", eval(e));
	}
	if (endl)
		print("\n");
	return N_EMPTY;
}

const Node* eval_if(span<const Node*> args) {
	if (args.size() != 2 && args.size() != 3)
		FAIL("(if) exactly 2 or 3 args");

	const Node* a = eval(args[0]);
	if (a == TRUE)
		return eval(args[1]);
	if (a == FALSE) {
		if (args.size() == 3)
			return eval(args[2]);
		return N_EMPTY;
	}
	if (is_error(a))
		return a;
	FAIL("(if) condition not bool", a);
}

const Node* eval_while(span<const Node*> args) {
	if (args.size() < 1)
		FAIL("(while) min 1 arg");

	auto cond = args[0];
	auto body = args.pop_front();
	while (true) {
		const Node* a = eval(cond);
		if (a == FALSE)
			return N_EMPTY;
		if (a != TRUE)
			FAIL("(while) condition not bool", a);
		// execute body
		for (auto e : body) {
			e = eval(e);
			if (e == N_BREAK)
				return N_EMPTY;
			if (e == N_CONTINUE)
				break;
			if (is_error(e))
				return e;
		}
	}
}

const Node* eval_break(span<const Node*> args) {
	if (args.size() != 0)
		FAIL("(break) no args");
	return N_BREAK;
}

const Node* eval_continue(span<const Node*> args) {
	if (args.size() != 0)
		FAIL("(continue) no args");
	return N_CONTINUE;
}

const Node* eval_equal(span<const Node*> args) {
	if (args.size() != 2)
		FAIL("(=) exactly 2 args");

	const Node* a = args[0];
	if (!a->is_strict_symbol())
		FAIL("(=) lhs must be symbol", a);
	const Node* b = eval(args[1]);
	if (is_error(b))
		return b;
	g_variables[a] = b;
	return b;
}

const Node* eval_def(const Node* n) {
	auto args = n->args();
	if (args.size() < 2)
		FAIL("(def) min 2 args");

	const Node* name = args[0];
	if (!name->is_strict_symbol())
		FAIL("(def) function name must be symbol", name);
	if (g_variables.count(name) != 0)
		FAIL("(def) function already defined", name);

	const Node* f_args = args[1];
	if (!f_args->is_list())
		FAIL("(def) function args must be list of symbols", f_args);
	for (auto e : f_args->elements())
	if (!e->is_strict_symbol())
		FAIL("(def) function arg must be symbol", e);

	g_variables[name] = n; // TODO convert to lambda instead
	return n;
}

const Node* eval_invoke(const Node* def, span<const Node*> args) {
	print("eval_invoke!\n");
	// TODO check that number of args match the definition
	// install args as variables
	// execute function body
	// (check for early return)
	// restore variables that were overwritten by args
	// return function result
	return N_EMPTY;
}

static const Node* eval(const Node* n) {
	if (n->is_string() || n->is_integer() || n == TRUE || n == FALSE)
		return n;
	if (n->is_keyword())
		FAIL("can't eval keyword", n);
	if (n->is_symbol()) {
		auto it = g_variables.find(n);
		if (it != g_variables.end())
			return it->second;
		FAIL("undefined symbol", n);
	}
	ASSERT_ALWAYS(n->is_list());
	// TODO empty list should be a symbol instead to avoid this check
	if (n->size() == 0)
		FAIL("list must begin with symbol", n);
	auto op = n->first();
	// TODO make these symbols into constexpr to be able to use switch
	if (op == EQUAL)
		return eval_equal(n->args());
	if (op == NOT)
		return eval_not(n->args());
	if (op == AND)
		return eval_and(n->args());
	if (op == OR)
		return eval_or(n->args());

	if (op == IF)
		return eval_if(n->args());
	if (op == WHILE)
		return eval_while(n->args());
	if (op == BREAK)
		return eval_break(n->args());
	if (op == CONTINUE)
		return eval_continue(n->args());

	if (op == PLUS)
		return eval_plus(n->args());
	if (op == PLUS_EQUAL)
		return eval_plus_equal(n->args());
	if (op == MINUS)
		return eval_minus(n->args());
	if (op == MINUS_EQUAL)
		return eval_minus_equal(n->args());
	if (op == STAR)
		return eval_star(n->args());
	if (op == STAR_EQUAL)
		return eval_star_equal(n->args());
	if (op == DIV_DIV)
		return eval_div_div(n->args());
	if (op == DIV_DIV_EQUAL)
		return eval_div_div_equal(n->args());
	if (op == PERCENT)
		return eval_percent(n->args());
	if (op == PERCENT_EQUAL)
		return eval_percent_equal(n->args());

	if (op == LESS_LESS)
		return eval_less_less(n->args());
	if (op == GREATER_GREATER)
		return eval_greater_greater(n->args());

	if (op == LESS)
		return eval_less(n->args());
	if (op == LESS_EQUAL)
		return eval_less_equal(n->args());
	if (op == GREATER)
		return eval_greater(n->args());
	if (op == GREATER_EQUAL)
		return eval_greater_equal(n->args());
	if (op == EQUAL_EQUAL)
		return eval_equal_equal(n->args());
	if (op == NOT_EQUAL)
		return eval_not_equal(n->args());

	if (op == PRINT)
		return eval_print(n->args(), false);
	if (op == PRINTLN)
		return eval_print(n->args(), true);
	if (op == DEF)
		return eval_def(n);
	if (op == NUMBER)
		return N_EMPTY;

	if (!op->is_symbol())
		FAIL("list must begin with symbol", n);

	auto it = g_variables.find(op);
	if (it != g_variables.end()) {
		auto d = it->second;
		if (d->is_list() && d->first() == DEF) {
			return eval_invoke(d, n->args());
		}
		FAIL("can't invoke", d);
	}

	FAIL("invalid op", op);
}

static const Node* eval_block(const Node* n) {
	if (!n->is_list())
		FAIL("list expected", n);
	const Node* result = N_EMPTY;
	for (auto e : n->elements()) {
		result = eval(e);
		if (is_error(result))
			break;
	}
	return result;
}

// TODO escape sequences in strings
// TODO in operator for maps
// TODO in operator for lists
// TODO in operator for strings
// TODO literal for lists [value value]
// TODO literal for maps {key value key value}
// TODO bigint mode for integers
// TODO big rational type
// TODO big decimal type
// TODO hexadecimal literals
// TODO binary literals
// TODO unit tests
// TODO bit operators
// TODO to_str, to_int, to_real
// TODO double datatype
// TODO function call
// TODO lambda
// TODO return from function
// TODO read terminal input
// TODO read file into string
// TODO write string into file
// TODO constants
// TODO list as data type
// TODO size of list / string
// TODO random access operator on list
// TODO string as list
// TODO operator + += on list
// TODO push_back, pop_back, pop_front, push_front
// TODO map (struct is just a map of fields and values)
// TODO break and continue without parenthesises
// TODO type operator (return type of variable / constant)
// TODO bootstrap!
// TODO for each loop
// TODO memory leaks?
// TODO importing public functions from other files

/*
int main(int argc, char** argv) {
	InitSegvHandler();

string_view code = R"END(
"This is comment"
(# This is comment)
(def sum (a b)
	(+ a b)
)
(set a "banana" 0)
(println "a[banana] = " (get a "banana"))
(println "(sum 4 5) -> " (sum 4 5))
(= names (list 1 3 5 7 9))
(println "names[2] = " (get names 2))
(for e names
	(println e)
)
(println (and true false))
(println (and true true))
(= a (+ 100 20))
(println a)
(println "hello world" "!")
(println (+ (* 3 4) -7))
(println (== 0 (+ -4 4)))
(if (< 3 2) (println "more") (println "less"))
(+ 1 4)
(println (- 2))
(= a 1)
(= b 2)
(while true
	(if (> a 10000) (break))
	(println a)
	(= c (+ a b))
	(= a b)
	(= b c)
)
)END";
	const Node* n = parse(code);
	if (n == nullptr)
		return 1;
	print("code: %s\n", n);
	const Node* e = eval_block(n);
	print("eval: %s\n", e);
	return 0;
}*/
