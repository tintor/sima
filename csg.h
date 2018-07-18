#pragma once
#include "triangle.h"

imesh3 csg_and(const imesh3& a, const imesh3& b);
imesh3 csg_or(const imesh3& a, const imesh3& b);
imesh3 csg_diff(const imesh3& a, const imesh3& b);
imesh3 csg_xor(const imesh3& a, const imesh3& b);
//imesh3 csg_cut(const imesh3& a, iplane3 p);
