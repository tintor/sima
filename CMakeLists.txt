cmake_minimum_required(VERSION 3.5)
project(sima)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
	
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include(FetchContent)

add_definitions(-O3 -mavx2 -mfma -ffast-math -ffp-contract=fast -DNDEBUG -Wfatal-errors -g)
#add_definitions(-O3 -mavx2 -lgmpxx -lgmp)
#add_definitions(-O3 -flto -DNDEBUG -mavx2 -lgmpxx -lgmp)

#  8.7s for (empty)
# 14.7s for -O3 -DNDEBUG -mavx2
# 16.1s for -O3 -DNDEBUG -mavx2 -flto

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ferror-limit=1")

# CUDA
# ----
if (NOT APPLE)
  find_package(CUDA 10.1 REQUIRED)
  add_definitions(-DHAVE_CUDA)
  add_definitions(--cuda-gpu-arch=sm_75 -Wno-unused-command-line-argument)
  include_directories(/usr/local/cuda-10.1/targets/x86_64-linux/include)
  link_directories(/usr/local/cuda-10.1/targets/x86_64-linux/lib)
  set(LIBRARIES ${LIBRARIES} ${CUDA_LIBRARIES} cublas cudnn cudart_static)
endif()

# Intel IPP
# ---------
add_definitions(-DHAVE_IPP)
include_directories(/opt/intel/ipp/include)
if(APPLE)
  link_directories(/opt/intel/ipp/lib)
else()
  link_directories(/opt/intel/ipp/lib/intel64)
endif()
set(LIBRARIES ${LIBRARIES} ippi ipps ippvm ippcore)

# Abseil
# ------
FetchContent_Declare(
    abseil-cpp
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp
    GIT_TAG 20200225.2
)
FetchContent_MakeAvailable(abseil-cpp)
include_directories(_deps/abseil-cpp-src)
set(LIBRARIES ${LIBRARIES}
absl::algorithm
absl::base
absl::debugging
absl::flat_hash_map
absl::flags
absl::memory
absl::meta
absl::numeric
absl::random_random
absl::strings
absl::synchronization
absl::time
absl::utility)

# {fmt}
# -----
#FetchContent_Declare(
#    fmt
#    GIT_REPOSITORY https://github.com/fmtlib/fmt
#    GIT_TAG 6.0.0
#)
#FetchContent_MakeAvailable(fmt)

# std::ranges
# -----------
FetchContent_Declare(
    ranges
    GIT_REPOSITORY https://github.com/ericniebler/range-v3
    GIT_TAG 0.10.0
)
FetchContent_MakeAvailable(ranges)
include_directories(_deps/ranges-src/include)

# ArrayFire
# ---------
# arrayfire.org/download
find_package(ArrayFire REQUIRED)
set(LIBRARIES ${LIBRARIES} ArrayFire::af)
include_directories(/opt/arrayfire/include)

# Klein - 3d geometric algebra
# ----------------------------
FetchContent_Declare(
    klein
    GIT_REPOSITORY https://github.com/jeremyong/Klein.git
    GIT_TAG v2.2.1
)
FetchContent_MakeAvailable(klein)
set(LIBRARIES ${LIBRARIES} klein::klein_sse42)

# GLM
# ---
# Apple: brew install glm
# Linux: sudo apt install libglm-dev

# Freetype
# --------
# Apple: brew install freetype
# Linux: sudo apt install libfreetype6-dev
find_package(Freetype REQUIRED)
if (APPLE)
	include_directories(/usr/local/Cellar/freetype/2.10.1/include/freetype2)
else()
	include_directories(/usr/include/freetype2)
endif()
set(LIBRARIES ${LIBRARIES} ${FREETYPE_LIBRARIES})

# glad
# ----
if (NOT APPLE)
	include_directories(glad/include)
endif()

# GLFW
# ----
# Apple: brew install glfw3
# Linux: sudo apt install libglfw3-dev libglfw3
find_package(glfw3 3.3 REQUIRED)
set(LIBRARIES ${LIBRARIES} glfw ${GLFW_LIBRARIES})
if (APPLE)
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework OpenGl")
endif()

# pthread and dl
# --------------
if (NOT APPLE)
	set(LIBRARIES ${LIBRARIES} dl pthread)
endif()

# lodepng
# -------
include_directories(lodepng)

# Eigen
# -----
# Apple: brew install eigen
# Linux: sudo apt install libeigen3-dev
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if (APPLE)
	include_directories(/usr/local/Cellar/eigen/3.3.7/include/eigen3)
else()
	include_directories(/usr/include/eigen3)
endif()
set(LIBRARIES ${LIBRARIES} Eigen3::Eigen)

# Boost
# -----
# Apple: brew install boost
# Linux: sudo apt install libboost-all-dev
find_package(Boost 1.67 REQUIRED)
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# ---------------------

add_subdirectory(core)
add_subdirectory(geom)
add_subdirectory(view)
add_subdirectory(sokoban)
add_subdirectory(santorini)

add_executable(test
	test.cc
	lisp.cc
	test_lisp.cc
	test_matrix_multiply.cc
	sim/test_integration.cc
	sim/test_gauss_seidel.cc)
target_link_libraries(test core ${LIBRARIES})

add_executable(arena arena.cc)
add_executable(main main.cc)
add_executable(piano_mover piano_mover.cc)
add_executable(collision collision.cc)
add_executable(run_gravity run_gravity.cc)
add_executable(pid pid.cc)

target_link_libraries(arena       ${LIBRARIES} core view)
target_link_libraries(main        ${LIBRARIES} core geom view)
target_link_libraries(piano_mover ${LIBRARIES} core geom view)
target_link_libraries(collision   ${LIBRARIES} core geom view)
target_link_libraries(run_gravity ${LIBRARIES} core geom view)
target_link_libraries(pid         ${LIBRARIES} core view)
