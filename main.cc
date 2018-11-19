#include "callstack.h"
#include "font.h"
#include "window.h"

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

int main(int argc, char** argv) {
	InitSegvHandler();

	auto window = CreateWindow({.width=800, .height=600, .title="Sima"});
	if (!window)
		return -1;

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);

	fprintf(stderr, "OpenGL version: [%s]\n", glGetString(GL_VERSION));
    glViewport(0, 0, 800, 600);
	glClearColor(0.0, 0.5, 0.0, 1.0);

	glEnable(GL_CULL_FACE);

	Font fontTNR("Times New Roman");
	Font fontArial("Arial");
	Font fontMonaco("/System/Library/Fonts/Monaco.dfont");
	if (!fontTNR.ok() || !fontArial.ok() || !fontMonaco.ok()) {
		return -1;
	}

	double x = 0;
	RunEventLoop(window, [&]() {
		x += 0.03;
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

    	glEnable(GL_BLEND);
    	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		fontTNR.render("This is Times New Roman text", x, 25, 1, "7FE030");
        fontArial.render("This is Arial text", 540, 570, 0.5, "40B0E0");
		fontMonaco.render("This is Monaco text", x, 100, 1, "FFFF00");
		glDisable(GL_BLEND);
	});
	return 0;
}
