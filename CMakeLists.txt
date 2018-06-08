#    LDFLAGS:  -L/usr/local/opt/llvm/lib
#    CPPFLAGS: -I/usr/local/opt/llvm/include

cmake_minimum_required(VERSION 2.8)

project(sima)

add_executable(sima main.cc rendering.cc lodepng/lodepng.cc)
add_executable(test primitives_test.cc primitives.cc integration_test.cc)

add_definitions(-g -O3 -flto -DNDEBUG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -ferror-limit=1")

add_subdirectory(glfw)
include_directories(glfw/include)
include_directories(lodepng)

if (APPLE)
    target_link_libraries(sima glfw ${GLFW_LIBRARIES})
else()
    find_package(GLEW REQUIRED)
    target_link_libraries(sima glfw dl ${GLEW_LIBRARIES} ${GLFW_LIBRARIES})
endif()
