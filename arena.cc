#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>

#include <sim/integration.h>
#include <sim/system.h>

#include <geom/classify.h>
#include <geom/triangle.h>
#include <geom/properties.h>
#include <geom/pose.h>

#include <view/glm.h>
#include <view/font.h>
#include <view/window.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>

bool gSimulate = false;
bool gSimulateTick = false;
bool gGravity = true;
bool gAirDrag = false;
bool gFriction = false;


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	const char* key_name = glfwGetKeyName(key, 0);
	print("key_callback [%s] key:%s scancode:%s action:%s mods:%s\n", key_name, key, scancode, action, mods);

	if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
		glfwSetWindowShouldClose(window, GL_TRUE);
		return;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_SPACE && mods == 0) {
		gSimulate ^= 1;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_SPACE && mods == GLFW_MOD_SHIFT) {
		gSimulateTick ^= 1;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_G && mods == 0) {
		gGravity ^= 1;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_A && mods == 0) {
		gAirDrag ^= 1;
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

constexpr int Width = 1200, Height = 900;

int main(int argc, char** argv) {
	InitSegvHandler();

	auto window = CreateWindow({.width=Width, .height=Height, .resizeable=false});
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, scroll_callback);

	glClearColor(0.0, 0.5, 0.0, 1.0);

	Shader shader(R"END(
		#version 330 core
		layout (location = 0) in vec2 pos;
		uniform mat3 transform;

		void main() {
		    vec3 p = transform * vec3(pos, 1.0);
		    gl_Position = vec4(p.x, p.y, 0.0, 1.0);
		}

		#version 330 core
		out vec4 color;

		void main() {
		    color = vec4(1.0, 0.5, 0.2, 1.0);
		}
	)END");
	UNIFORM(mat3, transform);

	VertexBuffer_vec2 buffer(25);

	Font timesNewRoman("Times New Roman");
	Font arial("Arial");
	Font monaco("/System/Library/Fonts/Monaco.dfont");

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

	RunEventLoop(window, [&]() {
		glClear(GL_COLOR_BUFFER_BIT);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		/*for (Body& body : bodies) {
			vec2 p = body.pos;
			p.x = p.x / Width * 800;
			p.y = p.y / Height * 600;
			monaco.render(
				format("pos (%f)\nvec (%f)", body.pos, body.vel),
				p.x, p.y, 0.3, "7FE030");
		}
		monaco.render(
			format("energy %s, spread %s, time %s", energy, energyMax - energyMin, time),
			5, 5, 0.3, "7FE030");*/
		glDisable(GL_BLEND);
	});
	return 0;
}
