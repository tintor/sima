#pragma once
#include "triangle.h"

mesh3 load_stl(string_view filename);
mesh3 load_ply(string_view filename);
