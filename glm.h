#pragma once

#define GLM_FORCE_RADIANS
#define GLM_FORCE_AVX2
#define GLM_ENABLE_EXPERIMENTAL
#define GLM_FORCE_CXX14
#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtx/quaternion.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_access.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/gtx/transform.hpp"

#include "common.h"
#include "format.h"

using glm::ivec2;
using glm::ivec3;
using lvec3 = glm::tvec3<long>;
using llvec3 = glm::tvec3<int128>;

using glm::vec3;
using glm::vec4;
using glm::dvec3;
using glm::dvec4;

using glm::mat3;
using glm::mat4;
using glm::dmat3;
using glm::dmat4;

using glm::quat;
using glm::dquat;

namespace glm {

template<typename T>
inline void format_e(std::string& s, std::string_view spec, glm::tvec2<T> v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
}

template<typename T>
inline void format_e(std::string& s, std::string_view spec, glm::tvec3<T> v) {
	::format_e(s, "", v.x);
	s += ' ';
	::format_e(s, "", v.y);
	s += ' ';
	::format_e(s, "", v.z);
}

}

namespace std {

template<typename T>
struct hash<glm::tvec2<T>> {
	size_t operator()(const glm::tvec2<T>& v) const {
		size_t seed = 0;
		hash_combine(seed, v.x);
		hash_combine(seed, v.y);
		return seed;
	}
};

template<typename T>
struct hash<glm::tvec3<T>> {
	size_t operator()(const glm::tvec3<T>& v) const {
		size_t seed = 0;
		hash_combine(seed, v.x);
		hash_combine(seed, v.y);
		hash_combine(seed, v.z);
		return seed;
	}
};

}
