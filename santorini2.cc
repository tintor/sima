#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>

#include <view/glm.h>
#include <view/font.h>
#include <view/window.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	const char* key_name = glfwGetKeyName(key, 0);
	print("key_callback [%s] key:%s scancode:%s action:%s mods:%s\n", key_name, key, scancode, action, mods);

	if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
		glfwSetWindowShouldClose(window, GL_TRUE);
		return;
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

constexpr int Width = 1000, Height = 1000;

int main(int argc, char** argv) {
	InitSegvHandler();

	auto window = CreateWindow({.width=Width, .height=Height, .resizeable=false});
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, scroll_callback);

	glClearColor(0.0, 0.0, 0.0, 1.0);

	Shader shader(R"END(
		#version 330 core
		layout (location = 0) in vec2 pos;
		layout (location = 1) in vec4 rgba;
		out vec4 pixel_rgba;
		uniform mat3 transform;

		void main() {
		    vec3 p = transform * vec3(pos, 1.0);
		    gl_Position = vec4(p.x, p.y, 0.0, 1.0);
            pixel_rgba = rgba;
		}

		#version 330 core
		in vec4 pixel_rgba;
		out vec4 color;

		void main() {
		    color = pixel_rgba;
		}
	)END");
	UNIFORM(mat3, transform);

	Font mono("JetBrainsMono-Medium.ttf");

	mat3 ortho;
	ortho[0][0] = 2.0f / Width;
	ortho[1][1] = 2.0f / Height;
	ortho[2][2] = 1.0f;
	ortho[2][0] = -1.0f;
	ortho[2][1] = -1.0f;

	VertexBuffer_vec2_rgba buffer(25);

	RunEventLoop(window, [&]() {
		glClear(GL_COLOR_BUFFER_BIT);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		mono.render(format("Hello world,\nMarko!"), 100, 100, 0.2, "FFFFFF");
		glDisable(GL_BLEND);

		glUseProgram(shader);
		buffer.bind();
		mat3 transform(ortho);
		transformUniform = transform;

        buffer.add({50, 50}, 0xFFFFFFFF);
        buffer.add({100, 50}, 0xFFFFFFFF);
        buffer.add({100, 100}, 0xFFFFFFFF);
        buffer.add({50, 100}, 0xFFFFFFFF);
		buffer.draw(GL_LINE_LOOP);
	});
	return 0;
}
