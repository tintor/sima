add_library(view
	window.cc
	font.cc
	rendering.cc
	shader.cc
	../lodepng/lodepng.cc)

if (NOT APPLE)
    add_library(view ../glad/src/glad.c)
endif()

