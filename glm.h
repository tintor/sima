#pragma once

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_transform_2d.hpp>
#include <glm/gtc/type_ptr.hpp>

using glm::dvec4;
using glm::dvec3;
using glm::dvec2;
using glm::vec2;
using glm::mat3;

namespace glm {
inline void format_e(string& s, string_view spec, dvec2 v) {
	::format_e(s, spec, v.x);
	s += ' ';
	::format_e(s, spec, v.y);
}

inline void format_e(string& s, string_view spec, dvec3 v) {
	::format_e(s, spec, v.x);
	s += ' ';
	::format_e(s, spec, v.y);
	s += ' ';
	::format_e(s, spec, v.z);
}
}
