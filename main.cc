// draw polygons with mouse
// render velocity and ang_velocity vectors
// move and rotate objects with mouse
// tangental friction
// notification messages that disappear
// convex shapes instead of balls
// concave shapes instead of convex
// hinges and 3 dof arm (turn off collision detection between objects that have hinge)
// joint friction
// PIC controller for the arm
// click on stop to move tip of end-effector there
// WASD to move end-effector
// control for each joint individually
// open/close end-effector

#include "callstack.h"
#include "font.h"
#include "window.h"
#include "shader.h"
#include "format.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_transform_2d.hpp>
#include <glm/gtc/type_ptr.hpp>

using glm::dvec2;
using glm::vec2;
using glm::mat3;

bool gSimulate = false;
bool gGravity = true;
bool gAirDrag = false;
bool gFriction = false;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
		glfwSetWindowShouldClose(window, GL_TRUE);
		return;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_SPACE && mods == 0) {
		gSimulate ^= 1;
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

class VertexBuffer_vec2 {
public:
	VertexBuffer_vec2(uint count) {
		glGenVertexArrays(1, &m_vao);
		glGenBuffers(1, &m_vbo);

		glBindVertexArray(m_vao);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vec2) * count, nullptr, GL_DYNAMIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, false, sizeof(vec2), 0);
	}

	void write(span<vec2> vertices) {
		glBindVertexArray(m_vao);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
		glBufferSubData(GL_ARRAY_BUFFER, 0, vertices.size() * sizeof(vec2), vertices.data());
	}

	void draw(uint type, span<vec2> vertices) {
		glBindVertexArray(m_vao);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
		glBufferSubData(GL_ARRAY_BUFFER, 0, vertices.size() * sizeof(vec2), vertices.data());
		glDrawArrays(type, 0, vertices.size());
	}

	void bind() {
		glBindVertexArray(m_vao);
	}

private:
	uint m_vbo, m_vao;
};

struct Body {
	double radius;
	double mass;
	dvec2 pos;
	dvec2 velocity = vec2(0, 0);
	double angle = 0;
	double ang_velocity = 0;
};

void ResolveCollision(Body& a, Body& b, dvec2 d) {
	double ua = glm::dot(d, a.velocity);
	double ub = glm::dot(d, b.velocity);
	auto ma = a.mass;
	auto mb = b.mass;
	double va = (ua * (ma - mb) + 2 * mb * ub) / (ma + mb);
	double vb = (ub * (mb - ma) + 2 * ma * ua) / (ma + mb);
	a.velocity += d * (va - ua);
	b.velocity += d * (vb - ub);
}

int main(int argc, char** argv) {
	InitSegvHandler();

	constexpr int Width = 1200, Height = 900;
	auto window = CreateWindow({.width=Width, .height=Height, .title="Sima", .resizeable=false});
	if (!window)
		return -1;

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, scroll_callback);

	fprintf(stderr, "OpenGL version: [%s]\n", glGetString(GL_VERSION));
	glClearColor(0.0, 0.5, 0.0, 1.0);

	//glEnable(GL_CULL_FACE);

	constexpr string_view VERT = R"END(
#version 330 core
layout (location = 0) in vec2 pos;
uniform mat3 transform;

void main() {
    vec3 p = transform * vec3(pos, 1.0);
    gl_Position = vec4(p.x, p.y, 0.0, 1.0);
}
)END";

	constexpr string_view FRAG = R"END(
#version 330 core
out vec4 color;

void main() {
    color = vec4(1.0, 0.5, 0.2, 1.0);
}
)END";

	Shader shader;
	if (!shader.load(VERT, FRAG)) {
		printf("shader failed\n");
		return -1;
	}

	const unsigned int transformLoc = glGetUniformLocation(shader, "transform");

	VertexBuffer_vec2 buffer(25);

	Font fontTNR("Times New Roman");
	Font fontArial("Arial");
	Font fontMonaco("/System/Library/Fonts/Monaco.dfont");
	if (!fontTNR.ok() || !fontArial.ok() || !fontMonaco.ok()) {
		return -1;
	}

	std::array<vec2, 25> v;
	for (int i = 0; i < v.size() - 1; i++) {
		auto a = 2 * M_PI / (v.size() - 1) * i;
		v[i].x = cos(a);
		v[i].y = sin(a);
	}
	v[24] = vec2(0, 0);
	buffer.write(v);

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

	std::vector<Body> bodies;
	bodies.resize(3);
	bodies[0] = {.mass=130*130, .radius=130, .pos=vec2(200, 200)};
	bodies[1] = {.mass=80*80, .radius=80, .pos=vec2(400, 400)};
	bodies[2] = {.mass=50*50, .radius=50, .pos=vec2(600, 600)};

	mat3 ortho;
	ortho[0][0] = 2.0f / Width;
	ortho[1][1] = 2.0f / Height;
	ortho[2][2] = 1.0f;
	ortho[2][0] = -1.0f;
	ortho[2][1] = -1.0f;

	dvec2 gravity(0, -10); // mm/s^2

	double time = 0;
	double energy = 0;
	for (Body& body : bodies) {
		energy += -glm::dot(body.pos, gravity) * body.mass + glm::dot(body.velocity, body.velocity) * body.mass / 2;
	}
	double energyMin = energy, energyMax = energy;

	// 250s -> 527 float vs 524 double
	RunEventLoop(window, [&]() {
		glClear(GL_COLOR_BUFFER_BIT);

		constexpr double dt = 0.01;
		if (gSimulate) {
			for (uint i = 0; i < bodies.size(); i++) {
				for (uint j = i + 1; j < bodies.size(); j++) {
					Body& a = bodies[i];
					Body& b = bodies[j];
					dvec2 d = a.pos - b.pos;
					float r = a.radius + b.radius;
					if (glm::dot(d, d) <= r * r) {
						ResolveCollision(a, b, glm::normalize(d));
					}
				}
			}
			for (Body& body : bodies) {
				if (gGravity && body.pos.y > body.radius) {
					body.velocity += gravity * dt;
				}
				if (gAirDrag) {
					body.velocity *= 0.999;
				}
				body.pos += body.velocity * dt;

				//double k = gFriction ? 0.95f : 1.0f;
				if (body.pos.x < body.radius && body.velocity.x < 0) {
					body.velocity.x = -body.velocity.x;// * k;
				}
				if (body.pos.x > Width - body.radius && body.velocity.x > 0) {
					body.velocity.x = -body.velocity.x;// * k;
				}
				if (body.pos.y < body.radius && body.velocity.y < 0) {
					body.velocity.y = -body.velocity.y;// * k;
				}
				if (body.pos.y > Height - body.radius && body.velocity.y > 0) {
					body.velocity.y = -body.velocity.y;// * k;
				}
			}
			time += dt;
		}

		double energy = 0;
		for (Body& body : bodies) {
			energy += -glm::dot(body.pos, gravity) * body.mass + glm::dot(body.velocity, body.velocity) * body.mass / 2;
		}
		if (energy < energyMin) energyMin = energy;
		if (energy > energyMax) energyMax = energy;

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		fontMonaco.render(format("Energy %.0f, spread %.0f, time %.2f, %.4f", energy, energyMax - energyMin, time, (energyMax - energyMin) / time), 25, 25, 0.3, "7FE030");
		glDisable(GL_BLEND);

		glUseProgram(shader);
		buffer.bind();
		for (const Body& body : bodies) {
			mat3 transform(ortho);
			transform = glm::translate(transform, vec2(body.pos));
			transform = glm::rotate(transform, float(body.angle));
			transform = glm::scale(transform, vec2(body.radius, body.radius));
			glUniformMatrix3fv(transformLoc, 1, false, glm::value_ptr(transform));
			glDrawArrays(GL_LINE_LOOP, 0, v.size());
		}
	});
	return 0;
}
