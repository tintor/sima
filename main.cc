#include "callstack.h"

// 3rd party
#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// std
#include <iostream>
using namespace std;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
		glfwSetWindowShouldClose(window, GL_TRUE);
		return;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_W && mods == 0) {
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_S && mods == 0) {
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_A && mods == 0) {
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_D && mods == 0) {
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
	if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
	}
}

void scroll_callback(GLFWwindow* window, double x, double y) {
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

void error_callback(int error, const char* message) {
	fprintf(stderr, "GLFW error %d: %s\n", error, message);
	exit(0);
}

void sigsegv_handler(int sig) {
	fprintf(stderr, "Error: signal %d:\n", sig);
	Callstack stack;
	string s;
	stack.write(s, {"sigsegv_handler(int)", "_sigtramp"});
	fputs(s.c_str(), stderr);
	exit(1);
}

int main(int argc, char** argv) {
	void sigsegv_handler(int sig);
	signal(SIGSEGV, sigsegv_handler);

    if (!glfwInit())
        return -1;

	glfwSetErrorCallback(error_callback);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(1000, 1000, "Sima", nullptr, nullptr);
	if (!window)
		return -1;
	glfwMakeContextCurrent(window);
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	    return -1;
	glfwSwapInterval(1); // VSYNC

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);

	fprintf(stderr, "OpenGL version: [%s]\n", glGetString(GL_VERSION));
	glViewport(0, 0, 1000, 1000);
	glClearColor(0.0, 0.5, 0.0, 1.0);

	// Blank screen workaround in OSX Mojave
	glfwPollEvents();
	int x, y;
	glfwGetWindowPos(window, &x, &y);
	glfwSetWindowPos(window, x+1, y);
	glfwSetWindowPos(window, x, y);

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);
		glDisable(GL_DEPTH_TEST);

		glfwSwapBuffers(window);
	}

	glfwTerminate();
	return 0;
}
