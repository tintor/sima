#    LDFLAGS:  -L/usr/local/opt/llvm/lib
#    CPPFLAGS: -I/usr/local/opt/llvm/include

cmake_minimum_required(VERSION 2.8)

project(sima)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(prod
	aabb.cc
	convex_hull.cc
	is_valid.cc
	tesselate.cc
	transform.cc
	properties.cc
	generators.cc
	solid_bsp_tree.cc
	csg.cc
	shape.cc
	file.cc
	primitives.cc
	timestamp.cc
	mesh_import.cc
	callstack.cc
	sphere.cc
	segment.cc
	classify.cc
	project.cc
	distance.cc
	broadphase.cc
	lodepng/lodepng.cc)

set(render
	font.cc
	rendering.cc
	glad/src/glad.c)

set(test
	test.cc
	test_span.cc
	test_aabb.cc
	test_util.cc
	test_aabb.cc
	test_classify.cc
	test_is_valid.cc
	test_transform.cc
	test_solid_bsp_tree.cc
	test_properties.cc
	test_generators.cc
	test_convex_hull.cc
	test_shape.cc
	test_primitives.cc
	test_tesselate.cc
	test_integration.cc
 	test_auto.cc
	#test_array_set.cc
	test_bits.cc
	test_dynamic_array.cc
	test_format.cc
	test_range.cc
	test_map.cc
	test_sphere.cc
	test_broadphase.cc
	test_zip.cc)

add_executable(a.sima main.cc ${prod} ${render})
add_executable(a.test ${prod} ${test})

add_definitions(-O3 -mavx2 -mfma -DNDEBUG)
#add_definitions(-O3 -mavx2 -lgmpxx -lgmp)
#add_definitions(-O3 -flto -DNDEBUG -mavx2 -lgmpxx -lgmp)

#  8.7s for (empty)
# 14.7s for -O3 -DNDEBUG -mavx2
# 16.1s for -O3 -DNDEBUG -mavx2 -flto

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -ferror-limit=1")

add_subdirectory(glm)
include_directories(glm)

# brew install freetype
include_directories(/usr/local/Cellar/freetype/2.9.1/include/freetype2)

include_directories(glad/include)

add_subdirectory(glfw)
include_directories(glfw/include)

include_directories(lodepng)
include_directories(../libraries)

find_package(Freetype REQUIRED)
if (APPLE)
	set(RENDER_LIBS glfw ${GLFW_LIBRARIES} ${FREETYPE_LIBRARIES})
	#target_link_libraries(a.test glfw ${GLFW_LIBRARIES})
else()
    find_package(GLEW REQUIRED)
	set(RENDER_LIBS glfw dl ${GLEW_LIBRARIES} ${GLFW_LIBRARIES} ${FREETYPE_LIBRARIES})
	#target_link_libraries(a.test glfw dl ${GLEW_LIBRARIES} ${GLFW_LIBRARIES})
endif()

target_link_libraries(a.sima ${RENDER_LIBS})
