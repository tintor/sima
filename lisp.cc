#include <core/auto.h>
#include <core/exception.h>
#include <core/format.h>
#include <core/range.h>
#include <lisp.h>

static bool is_integer(string_view s) {
    if (s.size() > 0 && (s[0] == '-' || s[0] == '+')) s = s.substr(1);
    if (s.empty()) return false;
    if (s[0] == '0') return s.size() == 1;
    for (char c : s)
        if (c < '0' || c > '9') return false;
    return true;
}

static long parse_long(string_view s) {
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
        if (s[0] == '+') i += 1;
        while (i < s.size()) {
            char c = s[i];
            a = a * 10 + c - '0';
            i += 1;
        }
    }
    return a;
}

struct Node;

unordered_map<string_view, const Node*> g_keywords;

struct Context {
    unordered_map<string_view, const Node*> symbols;
    unordered_map<const Node*, const Node*> variables;
};

static Context* g_context = nullptr;

bool equal(const Node*, const Node*);
size_t hash(const Node*);

enum class Type : uint8_t { List, Keyword, Symbol, String, Map };

struct MapHash {
    auto operator()(const Node* a) const { return hash(a); }
};

struct MapEqual {
    bool operator()(const Node* a, const Node* b) const { return equal(a, b); }
};

using Map = unordered_map<const Node*, const Node*, MapHash, MapEqual>;

struct Node {
   private:
    Type _type;
    uint _size = 0;
    void* _data = nullptr;  // TODO allocate memory inline

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

   public:
    static const Node* keyword(string_view a) {
        Node* n = alloc_raw(a.size());
        n->_type = Type::Keyword;
        memcpy(n->_data, a.data(), a.size());
        g_keywords.insert({n->symbol(), n});
        return n;
    }

    static const Node* symbol(string_view a) {
        auto it = g_context->symbols.find(a);
        if (it != g_context->symbols.end()) return it->second;

        Node* n = alloc_raw(a.size());
        n->_type = Type::Symbol;
        memcpy(n->_data, a.data(), a.size());
        g_context->symbols.insert({n->symbol(), n});
        return n;
    }

    static const Node* integer(long a);

    static const Node* str(string_view a) {
        Node* n = alloc_raw(a.size());
        n->_type = Type::String;
        memcpy(n->_data, a.data(), a.size());
        return n;
    }

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
        // TODO memcpy
        for (uint i : range(elements.size())) n->set(i, elements[i]);
        return n;
    }

    static const Node* list(std::initializer_list<const Node*> elements) {
        auto n = alloc_list(elements.size());
        // TODO memcpy
        for (uint i : range(elements.size())) n->set(i, elements.begin()[i]);
        return n;
    }

    auto type() const { return _type; }

    bool is_symbol() const { return _type == Type::Symbol; }
    bool is_keyword() const { return _type == Type::Keyword; }
    string_view symbol() const { return string_view((const char*)_data, _size); }

    bool is_integer() const { return reinterpret_cast<long>(this) & 1; }
    long integer() const { return reinterpret_cast<long>(this) >> 1; }

    bool is_string() const { return _type == Type::String; }
    string_view str() const { return string_view((const char*)_data, _size); }

    bool is_list() const { return _type == Type::List; }
    uint size() const { return _size; }
    const Node* first() const { return *(const Node**)_data; }
    auto elements() const { return span<const Node*>((const Node**)_data, _size); }
    auto args() const {
        const Node** data = (const Node**)_data;
        return span<const Node*>(data + 1, _size - 1);
    }
    const Node* get(uint index) const { return ((const Node**)_data)[index]; }
    void set(uint index, const Node* n) {
        const Node** data = (const Node**)_data;
        data[index] = n;
    }

    // Map
    static const Node* make_map() {
        Node* n = alloc_raw(sizeof(Map));
        n->_type = Type::Map;
        new (n->_data) Map();
        return n;
    }
    bool is_map() const { return _type == Type::Map; }
    const Map* map() const { return reinterpret_cast<const Map*>(_data); }
};

#define K(N, S) const Node* N = Node::keyword(#S);
K(NUMBER, #)
K(ERROR, error)
K(DEF, def)

K(TRUE, true)
K(FALSE, false)
K(AND, and)
K(OR, or)
K(NOT, not)

K(AMP, &)
K(PIPE, |)
K(CAPPA, ^)
K(TILDA, ~)

K(EQUAL, =)
K(PLUS, +)
K(PLUS_EQUAL, +=)
K(MINUS, -)
K(MINUS_EQUAL, -=)
K(STAR, *)
K(STAR_EQUAL, *=)
K(DIV, /)
K(DIV_EQUAL, /=)
K(PERCENT, %)
K(PERCENT_EQUAL, %=)

K(LESS_LESS, <<)
K(GREATER_GREATER, >>)

K(IF, if)
K(WHILE, while)
K(BREAK, break)
K(CONTINUE, continue)

K(LESS, <)
K(LESS_EQUAL, <=)
K(GREATER, >)
K(GREATER_EQUAL, >=)
K(EQUAL_EQUAL, ==)
K(NOT_EQUAL, !=)
K(EQUAL_EQUAL_EQUAL, == =)

K(PRINT, print)
K(PRINTLN, println)
K(PARSE, parse)
K(EVAL, eval)
K(SIZE, size)
K(MAP, map)
K(IN, in)
K(GET, get)
K(SET, set)

const Node* DIV_DIV = Node::keyword("//");
const Node* DIV_DIV_EQUAL = Node::keyword("//=");
const Node* EMPTY = Node::keyword("()");
const Node* QUOTE = Node::keyword("'");
#undef K

const Node* N_BREAK = Node::list(ERROR, BREAK);
const Node* N_CONTINUE = Node::list(ERROR, CONTINUE);

#define FAIL(...) return Node::list(ERROR, ##__VA_ARGS__)

const Node* Node::integer(long a) {
    constexpr long inline_max = numeric_limits<long>::max() >> 1;
    constexpr long inline_min = numeric_limits<long>::min() >> 1;
    if (a < inline_min || a > inline_max) FAIL("integer overflow");
    long v = (a << 1) | 1;
    return reinterpret_cast<const Node*>(v);
}

static bool space(char c) { return c == ' ' || c == '\t' || c == '\n'; }

static const Node* token_to_node(string_view s) {
    if (is_integer(s))
        // TODO return error if s can't fit into inline long
        return Node::integer(parse_long(s));
    auto it = g_keywords.find(s);
    if (it != g_keywords.end()) return it->second;
    return Node::symbol(s);
}

// TODO line no and column for syntax error
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
            if (size.size() == 1) FAIL("syntax: unexpected ')'");
            uint s = size.back();
            size.pop_back();
            const Node* n = (s == 0) ? EMPTY : Node::list(span<const Node*>(stack.data() + stack.size() - s, s));
            stack.resize(stack.size() - s);

            stack.push_back(n);
            size.back() += 1;
            p += 1;
            continue;
        }
        if (c == '"') {
            p += 1;
            int s = p;
            while (p < code.size() && code[p] != '"') p += 1;
            if (p == code.size()) FAIL("syntax: unterminated quote");
            const Node* n = Node::str(code.substr(s, p - s));
            stack.push_back(n);
            size.back() += 1;
            p += 1;
            continue;
        }
        // anything else is atom
        int s = p;
        while (p < code.size() && !space(code[p]) && code[p] != '(' && code[p] != ')') p += 1;
        stack.push_back(token_to_node(code.substr(s, p - s)));
        size.back() += 1;
    }
    if (size.size() != 1) FAIL("syntax: unterminated list");
    if (stack.size() == 0) FAIL("syntax: empty");
    if (stack.size() > 1) FAIL("syntax: only one top level expression expected");
    return stack[0];
}

int length(const Node* n) {
    ASSERT_ALWAYS(n != nullptr);
    if (n->is_integer()) {
        int s = 0;
        long c = n->integer();
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
            if (c == '"' || c == '\\') s += 1;
        return s + 2 + n->str().size();
    }
    if (n->is_symbol() || n->is_keyword()) return n->symbol().size();
    if (n->is_map()) {
        int s = 0;
        for (const auto& p : *n->map()) s += length(p.first) + length(p.second);
        return s + 5 + n->map()->size() * 2;
    }

    ASSERT_ALWAYS(n->is_list());
    if (n->size() == 0) return 2;
    int s = 0;
    for (const Node* e : n->elements()) s += length(e);
    return s + 1 + n->size();
}

void format_e(string& s, string_view spec, const Node* n);

void format_map(string& s, string_view spec, const Map& m) {
    s += "(map";
    for (const auto& p : m) {
        s += ' ';
        format_e(s, spec, p.first);
        s += ' ';
        format_e(s, spec, p.second);
    }
    s += ')';
}

void format_list(string& s, string_view spec, const Node* n) {
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
    for (int i : range(s_depth * 2)) s += ' ';
    s += "(\n";
    s_depth += 1;

    for (const Node* e : n->elements()) {
        ASSERT_ALWAYS(e != nullptr);
        if (!e->is_list() || length(e) <= 120) {
            for (int i : range(s_depth * 2)) s += ' ';
            format_e(s, spec, e);
        } else {
            format_e(s, spec, e);
        }
        s += '\n';
    }

    s_depth -= 1;
    for (int i : range(s_depth * 2)) s += ' ';
    s += ")";
}

void format_e(string& s, string_view spec, const Node* n) {
    ASSERT_ALWAYS(n != nullptr);
    if (n->is_integer()) {
        format_e(s, "%s", n->integer());
        return;
    }
    switch (n->type()) {
        case Type::Symbol:
        case Type::Keyword:
            s += n->symbol();
            return;
        case Type::String:
            s += '"';
            for (char c : n->str()) {
                if (c == '"' || c == '\\') s += '\\';
                s += c;
            }
            s += '"';
            return;
        case Type::Map:
            format_map(s, spec, *n->map());
            return;
        case Type::List:
            format_list(s, spec, n);
            return;
    }
}

static bool is_error(const Node* e) { return e->is_list() && e->size() > 0 && e->elements()[0] == ERROR; }

static const Node* eval(const Node* n);

const Node* eval_plus(span<const Node*> args) {
    if (args.size() == 2) {
        const Node* a = eval(args[0]);
        const Node* b = eval(args[1]);
        if (a->is_integer() && b->is_integer()) {
            long av = a->integer();
            long bv = b->integer();
            if (add_overflow(av, bv)) return Node::integer(cent(av) + cent(bv));
            return Node::integer(av + bv);
        }
        if (!a->is_integer() && is_error(a)) return a;
        if (!b->is_integer() && is_error(b)) return b;
    }
    if (args.size() == 0) FAIL("(+) min 1 arg");
    if (args.size() == 1) return eval(args[0]);

    cent acc = 0;
    for (const Node* e : args) {
        e = eval(e);
        if (!e->is_integer()) {
            if (is_error(e)) return e;
            FAIL("(+) arg not int", e);
        }
        acc += e->integer();
    }
    return Node::integer(acc);
}

const Node* eval_plus_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(+=) exactly 2 args");

    const Node* a = args[0];
    if (a->is_integer() || !a->is_symbol()) FAIL("(+=) lhs must be symbol", a);
    auto it = g_context->variables.find(a);
    if (it == g_context->variables.end()) FAIL("(+=) symbol not defined", a);
    auto& av = it->second;
    const Node* b = eval(args[1]);
    if (!b->is_integer()) {
        if (is_error(b)) return b;
        FAIL("(+=) arg not int", b);
    }
    b = Node::integer(av->integer() + b->integer());
    av = b;
    return b;
}

const Node* eval_minus(span<const Node*> args) {
    if (args.size() == 1) {
        auto a = eval(args[0]);
        if (!a->is_integer()) {
            if (is_error(a)) return a;
            FAIL("(-) arg not int", a);
        }
        return Node::integer(-a->integer());
    }

    if (args.size() == 2) {
        auto a = eval(args[0]);
        if (!a->is_integer()) {
            if (is_error(a)) return a;
            FAIL("(-) arg not int", a);
        }

        auto b = eval(args[1]);
        if (!b->is_integer()) {
            if (is_error(b)) return b;
            FAIL("(-) arg not int", b);
        }

        return Node::integer(a->integer() - b->integer());
    }
    FAIL("(-) exactly 1 or 2 args");
}

const Node* eval_minus_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(-=) exactly 2 args");

    const Node* a = args[0];
    if (!a->is_symbol()) FAIL("(-=) lhs must be symbol", a);
    auto it = g_context->variables.find(a);
    if (it == g_context->variables.end()) FAIL("(+=) symbol not defined", a);
    auto& av = it->second;
    const Node* b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(-=) arg not int", b);
    b = Node::integer(av->integer() - b->integer());
    av = b;
    return b;
}

const Node* eval_not(span<const Node*> args) {
    if (args.size() != 1) FAIL("(not) exactly 1 arg");

    auto a = eval(args[0]);
    if (!a->is_integer() && is_error(a)) return a;
    if (a == TRUE) return FALSE;
    if (a == FALSE) return TRUE;
    FAIL("(not) arg not bool", a);
}

const Node* eval_and(span<const Node*> args) {
    if (args.size() == 0) FAIL("(and) min 1 arg");

    for (auto e : args) {
        e = eval(e);
        if (!e->is_integer() && is_error(e)) return e;
        if (e == FALSE) return FALSE;
        if (e != TRUE) FAIL("(and) arg not bool", e);
    }
    return TRUE;
}

const Node* eval_or(span<const Node*> args) {
    if (args.size() == 0) FAIL("(or) min 1 arg");

    for (auto e : args) {
        e = eval(e);
        if (!e->is_integer() && is_error(e)) return e;
        if (e == TRUE) return TRUE;
        if (e != FALSE) FAIL("(or) arg not bool", e);
    }
    return FALSE;
}

const Node* eval_star(span<const Node*> args) {
    if (args.size() < 2) FAIL("(*) min 2 args");
    cent acc = 1;
    for (auto e : args) {
        e = eval(e);
        if (is_error(e)) return e;
        if (!e->is_integer()) FAIL("(*) arg not int", e);
        acc *= e->integer();
    }
    return Node::integer(acc);
}

const Node* eval_star_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(*=) exactly 2 args");

    const Node* a = args[0];
    if (!a->is_symbol()) FAIL("(*=) lhs must be symbol", a);
    auto it = g_context->variables.find(a);
    if (it == g_context->variables.end()) FAIL("(*=) symbol not defined", a);
    auto& av = it->second;
    const Node* b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(*=) arg not int", b);
    b = Node::integer(av->integer() * b->integer());
    av = b;
    return b;
}

const Node* eval_div_div(span<const Node*> args) {
    if (args.size() != 2) FAIL("(//) exactly 2 args");

    auto a = eval(args[0]);
    if (is_error(a)) return a;
    if (!a->is_integer()) FAIL("(//) arg not int", a);

    auto b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(//) arg not int", b);

    return Node::integer(a->integer() / b->integer());
}

const Node* eval_div_div_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(//=) exactly 2 args");

    const Node* a = args[0];
    if (!a->is_symbol()) FAIL("(//=) lhs must be symbol", a);
    auto it = g_context->variables.find(a);
    if (it == g_context->variables.end()) FAIL("(//=) symbol not defined", a);
    auto& av = it->second;
    const Node* b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(//=) arg not int", b);
    b = Node::integer(av->integer() / b->integer());
    av = b;
    return b;
}

const Node* eval_percent(span<const Node*> args) {
    if (args.size() != 2) FAIL("(%) exactly 2 args");

    auto a = eval(args[0]);
    if (is_error(a)) return a;
    if (!a->is_integer()) FAIL("(%) arg not int", a);

    auto b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(%) arg not int", b);

    return Node::integer(a->integer() % b->integer());
}

const Node* eval_percent_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(%=) exactly 2 args");

    const Node* a = args[0];
    if (!a->is_symbol()) FAIL("(%=) lhs must be symbol", a);
    auto it = g_context->variables.find(a);
    if (it == g_context->variables.end()) FAIL("(%=) symbol not defined", a);
    auto& av = it->second;
    const Node* b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(%=) arg not int", b);
    b = Node::integer(av->integer() % b->integer());
    av = b;
    return b;
}

const Node* eval_less_less(span<const Node*> args) {
    if (args.size() != 2) FAIL("(<<) exactly 2 args");

    auto a = eval(args[0]);
    if (is_error(a)) return a;
    if (!a->is_integer()) FAIL("(<<) arg not int", a);

    auto b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(<<) arg not int", b);

    return Node::integer(a->integer() << b->integer());
}

const Node* eval_greater_greater(span<const Node*> args) {
    if (args.size() != 2) FAIL("(>>) exactly 2 args");

    auto a = eval(args[0]);
    if (is_error(a)) return a;
    if (!a->is_integer()) FAIL("(>>) arg not int", a);

    auto b = eval(args[1]);
    if (is_error(b)) return b;
    if (!b->is_integer()) FAIL("(>>) arg not int", b);

    return Node::integer(a->integer() >> b->integer());
}

#define COMPARE(OP)                                                                             \
    if (args.size() != 2) FAIL("(" #OP ") exactly 2 args");                                     \
    const Node* a = eval(args[0]);                                                              \
    const Node* b = eval(args[1]);                                                              \
    if (a->is_integer() && b->is_integer()) return a->integer() OP b->integer() ? TRUE : FALSE; \
    if (!a->is_integer()) {                                                                     \
        if (is_error(a)) return a;                                                              \
        FAIL("(" #OP ") arg not int", a);                                                       \
    }                                                                                           \
    if (!b->is_integer()) {                                                                     \
        if (is_error(b)) return b;                                                              \
        FAIL("(" #OP ") arg not int", b);                                                       \
    }                                                                                           \
    THROW(runtime_error);

const Node* eval_less(span<const Node*> args) { COMPARE(<) }
const Node* eval_less_equal(span<const Node*> args) { COMPARE(<=) }
const Node* eval_greater(span<const Node*> args) { COMPARE(>) }
const Node* eval_greater_equal(span<const Node*> args) { COMPARE(>=) }

static bool raw_equal(const Node* a, const Node* b);

size_t hash(const Node* a) {
    if (a->is_integer()) return a->integer();

    size_t h = 0;
    switch (a->type()) {
        case Type::String:
            return a->str().size();  // TODO horible!
        case Type::Symbol:
        case Type::Keyword:
            return reinterpret_cast<size_t>(a);
        case Type::Map:
            for (const auto& p : *a->map()) h += hash(p.first) + hash(p.second);
            return h;
        case Type::List:
            for (auto e : a->elements()) h += hash(e);
            return h;  // TODO horible!
    }
}

bool equal(const Node* a, const Node* b) {
    if (a == b) return true;
    if (a->is_integer() || b->is_integer()) return false;
    return raw_equal(a, b);
}

static bool maps_equal(const Map& a, const Map& b) {
    if (a.size() != b.size()) return false;
    for (const auto& p : a) {
        auto q = b.find(p.first);
        if (q == b.end() || !equal(p.second, q->second)) return false;
    }
    return true;
}

// assumes a and b are not integers
static bool raw_equal(const Node* a, const Node* b) {
    if (a->type() != b->type()) return false;

    switch (a->type()) {
        case Type::String:
            return a->str() == b->str();
        case Type::Symbol:
            return a == b;
        case Type::Keyword:
            return a == b;
        case Type::Map:
            return maps_equal(*a->map(), *b->map());
        case Type::List:
            if (a->size() != b->size()) return false;
            for (auto i : range(a->size()))
                if (!equal(a->get(i), b->get(i))) return false;
            return true;
    }
    THROW(runtime_error);
}

const Node* eval_equal_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(==) exactly 2 args");
    auto a = eval(args[0]);
    if (a->is_integer()) {
        auto b = eval(args[1]);
        if (b->is_integer()) return a == b ? TRUE : FALSE;
        if (is_error(b)) return b;
        return FALSE;
    } else {
        if (is_error(a)) return a;
        auto b = eval(args[1]);
        if (b->is_integer()) return FALSE;
        return raw_equal(a, b) ? TRUE : FALSE;
    }
}

const Node* eval_not_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(!=) exactly 2 args");
    auto a = eval(args[0]);
    if (a->is_integer()) {
        auto b = eval(args[1]);
        if (b->is_integer()) return a != b ? TRUE : FALSE;
        if (is_error(b)) return b;
        return FALSE;
    } else {
        if (is_error(a)) return a;
        auto b = eval(args[1]);
        if (b->is_integer()) return FALSE;
        return !raw_equal(a, b) ? TRUE : FALSE;
    }
}

const Node* eval_print(span<const Node*> args, bool endl) {
    for (auto e : args) {
        e = eval(e);
        if (e->is_integer()) {
            print("%s", e->integer());
            continue;
        }
        if (is_error(e)) return e;
        if (e->is_string())
            print("%s", e->str());
        else
            print("%s", e);
    }
    if (endl) print("\n");
    return EMPTY;
}

const Node* eval_if(span<const Node*> args) {
    if (args.size() != 2 && args.size() != 3) FAIL("(if) exactly 2 or 3 args");

    const Node* a = eval(args[0]);
    if (a == TRUE) return eval(args[1]);
    if (a == FALSE) {
        if (args.size() == 3) return eval(args[2]);
        return EMPTY;
    }
    if (!a->is_integer() && is_error(a)) return a;
    FAIL("(if) condition not bool", a);
}

// TODO return last expression from while's body
const Node* eval_while(span<const Node*> args) {
    if (args.size() < 1) FAIL("(while) min 1 arg");

    auto cond = args[0];
    auto body = args.pop_front();
    while (true) {
        const Node* a = eval(cond);
        if (a == FALSE) return EMPTY;
        if (a != TRUE) {
            if (!a->is_integer() && is_error(a)) return a;
            FAIL("(while) condition not bool", a);
        }
        // execute body
        for (auto e : body) {
            e = eval(e);
            if (e == N_BREAK) return EMPTY;
            if (e == N_CONTINUE) break;
            if (!e->is_integer() && is_error(e)) return e;
        }
    }
}

const Node* eval_equal(span<const Node*> args) {
    if (args.size() != 2) FAIL("(=) exactly 2 args");

    const Node* a = args[0];
    if (a->is_integer() || !a->is_symbol()) FAIL("(=) lhs must be symbol", a);
    const Node* b = eval(args[1]);
    if (!b->is_integer() && is_error(b)) return b;
    g_context->variables[a] = b;
    return b;
}

const Node* eval_def(const Node* n) {
    auto args = n->args();
    if (args.size() < 2) FAIL("(def) min 2 args");

    const Node* name = args[0];
    if (!name->is_symbol()) FAIL("(def) function name must be symbol", name);
    if (g_context->variables.count(name) != 0) FAIL("(def) function already defined", name);

    const Node* f_args = args[1];
    if (!f_args->is_list()) FAIL("(def) function args must be list of symbols", f_args);
    for (auto e : f_args->elements())
        if (!e->is_symbol()) FAIL("(def) function arg must be symbol", e);

    g_context->variables[name] = n;  // TODO convert to lambda instead
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
    return EMPTY;
}

static const Node* eval(const Node* n) {
    if (n->is_integer() || n->is_string() || n == TRUE || n == FALSE || n == EMPTY) return n;

    if (n->is_keyword()) {
        if (n == BREAK) return N_BREAK;
        if (n == CONTINUE) return N_CONTINUE;
        FAIL("can't eval keyword", n);
    }

    if (n->is_symbol()) {
        auto it = g_context->variables.find(n);
        if (it != g_context->variables.end()) return it->second;
        FAIL("undefined symbol", n);
    }

    ASSERT_ALWAYS(n->is_list());
    auto op = n->first();
    if (op->is_integer() || op->is_string() || op == TRUE || op == FALSE || op == EMPTY) return n;

    // TODO make these symbols into constexpr to be able to use switch
    auto args = n->args();
    if (op == EQUAL) return eval_equal(args);
    if (op == NOT) return eval_not(args);
    if (op == AND) return eval_and(args);
    if (op == OR) return eval_or(args);

    if (op == IF) return eval_if(args);
    if (op == WHILE) return eval_while(n->args());

    if (op == PLUS) return eval_plus(n->args());
    if (op == PLUS_EQUAL) return eval_plus_equal(n->args());
    if (op == MINUS) return eval_minus(n->args());
    if (op == MINUS_EQUAL) return eval_minus_equal(n->args());
    if (op == STAR) return eval_star(n->args());
    if (op == STAR_EQUAL) return eval_star_equal(n->args());
    if (op == DIV_DIV) return eval_div_div(n->args());
    if (op == DIV_DIV_EQUAL) return eval_div_div_equal(n->args());
    if (op == PERCENT) return eval_percent(n->args());
    if (op == PERCENT_EQUAL) return eval_percent_equal(n->args());

    if (op == LESS_LESS) return eval_less_less(n->args());
    if (op == GREATER_GREATER) return eval_greater_greater(n->args());

    if (op == LESS) return eval_less(n->args());
    if (op == LESS_EQUAL) return eval_less_equal(n->args());
    if (op == GREATER) return eval_greater(n->args());
    if (op == GREATER_EQUAL) return eval_greater_equal(n->args());
    if (op == EQUAL_EQUAL) return eval_equal_equal(n->args());
    if (op == NOT_EQUAL) return eval_not_equal(n->args());

    if (op == PRINT) return eval_print(n->args(), false);
    if (op == PRINTLN) return eval_print(n->args(), true);
    if (op == DEF) return eval_def(n);
    if (op == NUMBER) return EMPTY;
    if (op == PARSE) {
        if (args.size() != 1) FAIL("(parse) exactly 1 arg");
        auto a = eval(args[0]);
        if (a->is_integer()) FAIL("(parse) arg not string", a);
        if (is_error(a)) return a;
        return parse(a->str());
    }
    if (op == EVAL) {
        if (args.size() != 1) FAIL("(eval) exactly 1 arg");
        auto a = eval(args[0]);
        if (a->is_integer() || is_error(a)) return a;
        return eval(a);
    }
    if (op == QUOTE) {
        if (args.size() == 0) return EMPTY;
        if (args.size() == 1) return args[0];
        return Node::list(args);
    }
    if (op == ERROR) return n;
    if (op->is_list()) {
        auto res = EMPTY;
        for (auto e : n->elements()) {
            auto a = eval(e);
            if (!a->is_integer() && is_error(a)) return a;
            res = a;
        }
        return res;
    }
    if (op == EQUAL_EQUAL_EQUAL) {
        if (args.size() != 2) FAIL("(===) exactly 2 args");
        auto a = eval(args[0]);
        if (!a->is_integer() && is_error(a)) return a;
        auto b = eval(args[1]);
        if (!b->is_integer() && is_error(b)) return b;
        return a == b ? TRUE : FALSE;
    }
    if (op == SIZE) {
        if (args.size() != 1) FAIL("(size) exactly 1 arg");
        auto a = eval(args[0]);
        if (a->is_integer()) FAIL("(size) arg not string, symbol or list", a);
        if (is_error(a)) return a;
        if (a->is_string()) return Node::integer(a->str().size());
        if (a->is_map()) return Node::integer(a->map()->size());
        if (a->is_symbol()) return Node::integer(a->symbol().size());
        if (a == EMPTY) return Node::integer(0);
        if (a->is_list()) return Node::integer(a->size());
        FAIL("(size) arg not string, symbol, list or map", a);
    }
    if (op == MAP) {
        if (args.size() % 2 != 0) FAIL("(map) even numbers of args");
        auto n = Node::make_map();
        auto& m = *const_cast<Map*>(n->map());
        for (auto i : range(args.size() / 2)) {
            auto k = eval(args[i * 2]);
            // key can be: list, string, symbol, integer, bool
            if (!k->is_integer() && is_error(k)) return k;
            auto v = eval(args[i * 2 + 1]);
            if (!v->is_integer() && is_error(v)) return v;
            m[k] = v;
        }
        return n;
    }
    if (op == IN) {
        if (args.size() != 2) FAIL("(in) exactly 2 args");
        auto a = eval(args[0]);
        if (!a->is_integer() && is_error(a)) return a;
        auto b = eval(args[1]);
        if (b->is_integer()) FAIL("(in) second arg must be map", b);
        if (is_error(b)) return b;
        if (!b->is_map()) FAIL("(in) second arg must be map", b);

        return b->map()->count(a) > 0 ? TRUE : FALSE;
    }
    if (op == GET) {
        if (args.size() != 2) FAIL("(get) exactly 2 args");

        auto a = eval(args[0]);
        if (a->is_integer()) FAIL("(get) first arg must be list or map", a);
        if (is_error(a)) return a;
        auto b = eval(args[1]);
        if (!b->is_integer() && is_error(b)) return b;

        if (a->is_list()) {
            if (!b->is_integer()) FAIL("(get) second arg must be integer for list", b);
            long i = b->integer();
            if (i < 0 || i >= a->size()) return EMPTY;
            return a->get(i);
        }
        if (a->is_map()) {
            auto it = a->map()->find(b);
            if (it != a->map()->end()) return it->second;
            return EMPTY;
        }
        FAIL("(get) first arg must be list or map", a);
    }
    if (op == SET) {
        if (args.size() != 3) FAIL("(set) exactly 3 args");

        auto a = eval(args[0]);
        if (a->is_integer()) FAIL("(set) first arg must be list or map", a);
        if (is_error(a)) return a;
        auto b = eval(args[1]);
        if (!b->is_integer() && is_error(b)) return b;
        auto c = eval(args[2]);
        if (!c->is_integer() && is_error(c)) return c;

        if (a->is_list()) {
            if (!b->is_integer()) FAIL("(set) second arg must be integer for list", b);
            long i = b->integer();
            if (i < 0 || i >= a->size()) FAIL("(set) index out of bounds", b);
            const_cast<Node*>(a)->set(i, c);
            return c;
        }
        if (a->is_map()) {
            auto& m = *const_cast<Map*>(a->map());
            m[b] = c;
            return c;
        }
        FAIL("(set) first arg must be list or map", a);
    }

    if (op->is_symbol()) {
        auto it = g_context->variables.find(op);
        if (it != g_context->variables.end()) {
            auto d = it->second;
            if (d->is_list() && d->first() == DEF) {
                return eval_invoke(d, n->args());
            }
            FAIL("can't invoke", d);
        }
    }

    FAIL("invalid op", op);
}

static const Node* eval_block(const Node* n) {
    if (!n->is_list()) FAIL("list expected", n);
    const Node* result = EMPTY;
    for (auto e : n->elements()) {
        result = eval(e);
        if (is_error(result)) break;
    }
    return result;
}

// TODO escape sequences in strings \" \\ \n \t
// TODO function call !!
// TODO lambda !!
// TODO closure !!
// TODO early return from function !!
// TODO recursion !!
// TODO macros (not eval-ing function args by default) !!
// TODO read terminal input
// TODO read file into string
// TODO write string into file
// TODO get char from string
// TODO resizeable lists! push_back element, pop_back element, append list to list
// TODO less parenthesis!
// TODO infix operators
// TODO importing public functions from other files -> start working on standard library
// TODO minimization: move stuff that can be moved to library
// TODO releasing memory! and catching leaks!

// TODO bootstrap!

// TODO tests for print
// TODO multiline strings
// TODO for each loop
// TODO type operator (return type of variable / constant)
// TODO pop_front, push_front
// TODO big integer type
// TODO big rational type
// TODO big decimal type
// TODO hexadecimal / binary literals
// TODO in operator for strings
// TODO literal for lists [value value]
// TODO literal for maps {key:value key:value}
// TODO bit operators
// TODO to_str, to_int, to_real
// TODO double datatype
// TODO complex double datatype
// TODO constants

Lisp::Lisp() { m_handle = new Context; }

Lisp::~Lisp() { delete reinterpret_cast<Context*>(m_handle); }

string Lisp::format(Lisp::node n) { return ::format("%s", reinterpret_cast<const Node*>(n)); }

Lisp::node Lisp::parse(string_view code) {
    g_context = reinterpret_cast<Context*>(m_handle);
    ON_SCOPE_EXIT(g_context = nullptr);
    auto n = ::parse(code);
    return reinterpret_cast<Lisp::node>(n);
}

Lisp::node Lisp::eval(Lisp::node n) {
    g_context = reinterpret_cast<Context*>(m_handle);
    ON_SCOPE_EXIT(g_context = nullptr);
    return reinterpret_cast<Lisp::node>(::eval(reinterpret_cast<const Node*>(n)));
}

bool Lisp::equal(Lisp::node a, Lisp::node b) {
    return ::equal(reinterpret_cast<const Node*>(a), reinterpret_cast<const Node*>(b));
}
