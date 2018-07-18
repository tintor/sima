#pragma once

#include "triangle.h"
#include <string_view>

imesh3 load_stl(std::string_view filename);
imesh3 load_ply(std::string_view filename);
