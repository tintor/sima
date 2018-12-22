#    LDFLAGS:  -L/usr/local/opt/llvm/lib
#    CPPFLAGS: -I/usr/local/opt/llvm/include

cmake_minimum_required(VERSION 2.8)

project(sima)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(prod
	core/callstack.cc
	core/file.cc
	geom/aabb.cc
	geom/classify.cc
	geom/convex_hull.cc
	geom/csg.cc
	geom/distance.cc
	geom/is_valid.cc
	geom/gjk.cc
	tesselate.cc
	transform.cc
	properties.cc
	generators.cc
	solid_bsp_tree.cc
	shape.cc
	primitives.cc
	timestamp.cc
	mesh_import.cc
	buffer.cc
	sphere.cc
	segment.cc
	project.cc
	broadphase.cc
	lodepng/lodepng.cc)

set(render
	view/window.cc
	view/font.cc
	view/rendering.cc
	view/shader.cc
	glad/src/glad.c)

set(test
	test.cc
	core/test_array_deque.cc
	core/test_array_set.cc
 	core/test_auto.cc
	core/test_bits.cc
	core/test_dynamic_array.cc
	core/test_format.cc
	core/test_map.cc
	core/test_range.cc
	core/test_span.cc
	core/test_util.cc
	core/test_zip.cc
	sim/test_integration.cc
	geom/test_aabb.cc
	geom/test_buffer.cc
	geom/test_classify.cc
	geom/test_convex_hull.cc
	geom/test_is_valid.cc
	test_vector.cc
	test_transform.cc
	test_solid_bsp_tree.cc
	test_properties.cc
	test_generators.cc
	test_shape.cc
	test_primitives.cc
	test_tesselate.cc
	test_sphere.cc
	test_broadphase.cc
	test_quaternion.cc)

add_executable(main main.cc ${prod} ${render})
add_executable(run_gravity run_gravity.cc ${prod} ${render})
add_executable(test ${prod} ${test})
add_executable(pid pid.cc ${prod} ${render})

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

include_directories(.)
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

target_link_libraries(main ${RENDER_LIBS})
target_link_libraries(run_gravity ${RENDER_LIBS})
target_link_libraries(pid ${RENDER_LIBS})
