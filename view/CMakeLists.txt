if (NOT APPLE)
    add_library(view
	window.cc
	font.cc
	rendering.cc
	shader.cc
	../lodepng/lodepng.cc
    ../glad/src/glad.c)
else()
    add_library(view
	window.cc
	font.cc
	rendering.cc
	shader.cc
	../lodepng/lodepng.cc)
endif()

