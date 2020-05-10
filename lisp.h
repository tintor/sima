#pragma once
#include <string>
#include <string_view>

class Lisp {
   public:
    using node = const int*;

    Lisp();
    ~Lisp();

    std::string format(node n);
    node parse(std::string_view code);
    node eval(node n);
    static bool equal(node, node);

   private:
    void* m_handle;
};
