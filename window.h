#pragma once

// 3rd party
#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <functional>

struct WindowDef {
	// width / height is size of render area (doesn't include window decorations)
	int width = 800;
	int height = 600;
	const char* title = "";
	bool fullscreen = false;
	bool vsync = true;
};

GLFWwindow* CreateWindow(WindowDef wd);
void RunEventLoop(GLFWwindow* window, const std::function<void()>& callback);
