#pragma once
#include "mesh.h"

mesh3 csg_and(const mesh3& a, const mesh3& b);
mesh3 csg_or(const mesh3& a, const mesh3& b);
mesh3 csg_diff(const mesh3& a, const mesh3& b);
mesh3 csg_xor(const mesh3& a, const mesh3& b);
//imesh3 csg_cut(const imesh3& a, iplane3 p);
