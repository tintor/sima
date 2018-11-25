// collision between concave shapes
// rotation!
// hinges and 3 dof arm (turn off collision detection between objects that have hinge)
// draw polygons with mouse
// render velocity and ang_velocity vectors
// move and rotate objects with mouse
// tangental friction
// notification messages that disappear
// joint friction
// PID controller for the arm
// click on stop to move tip of end-effector there
// WASD to move end-effector
// control for each joint individually
// open/close end-effector

// 3D:
// balls
// convex shapes (no rotation)
// rotation
// hinges
// concave shapes

// general
#include "callstack.h"
#include "format.h"
#include "integration.h"
#include "util.h"
#include "glm.h"
#include "classify.h"
#include "triangle.h"
#include "properties.h"

// rendering
#include "font.h"
#include "window.h"
#include "shader.h"
#include "vertex_buffer.h"

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

struct Body {
	polygon2 shape;
	std::vector<vec2> shapef;

	aabb2 box;
	double radius;

	double mass;
	dvec2 pos;
	dvec2 velocity = vec2(0, 0);
	double angle = 0;
	double ang_velocity = 0;
};

void Interact(Body& a, Body& b) {
	dvec2 d = a.pos - b.pos;
	float r = a.radius + b.radius;
	if (glm::dot(d, d) > r * r)
		return;
	d = glm::normalize(d);

	vector<Contact2> contacts;
	int c = Classify(a.shape, b.shape, &contacts);
	if (c > 0)
		return;

	double ua = glm::dot(d, a.velocity);
	double ub = glm::dot(d, b.velocity);
	auto ma = a.mass;
	auto mb = b.mass;
	double va = (ua * (ma - mb) + 2 * mb * ub) / (ma + mb);
	double vb = (ub * (mb - ma) + 2 * ma * ua) / (ma + mb);
	a.velocity += d * (va - ua);
	b.velocity += d * (vb - ub);
}

void SetShape(Body& body, const polygon2& poly) {
	double2 com = CenterOfMass(poly);
	body.shape = poly;
	for (double2& v : body.shape)
		v -= com;

	body.shapef.resize(body.shape.size());
	for (uint i = 0; i < body.shape.size(); i++)
		body.shapef[i] = vec2(body.shape[i].x, body.shape[i].y);

	double r2 = 0;
	for (double2 v : body.shape)
		maximize(r2, dot(v, v));
	body.radius = sqrt(r2);

	body.box = aabb2(body.shape);
}

int main(int argc, char** argv) {
	InitSegvHandler();

	constexpr int Width = 1200, Height = 900;
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

	std::array<vec2, 25> v;
	for (int i = 0; i < v.size() - 1; i++) {
		auto a = 2 * M_PI / (v.size() - 1) * i;
		v[i].x = cos(a);
		v[i].y = sin(a);
	}
	v[24] = vec2(0, 0);
	std::array<vec2, 25> w;
	buffer.write(v);

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

	std::vector<Body> bodies;
	bodies.push_back({.mass=130*130, .radius=130, .pos=vec2(200, 200)});
	bodies.push_back({.mass=80*80, .radius=80, .pos=vec2(400, 400)});
	bodies.push_back({.mass=50*50, .radius=50, .pos=vec2(600, 600)});

	polygon2 shape;
	shape.push_back(double2{0, 0});
	shape.push_back(double2{200, 0});
	shape.push_back(double2{200, -100});
	shape.push_back(double2{100, -100});
	shape.push_back(double2{100, -200});
	shape.push_back(double2{0, -200});
	SetShape(bodies[0], shape);

	shape.clear();
	for (int i = 0; i < 24; i++) {
		auto a = 2 * M_PI / 24 * i;
		shape.push_back(double2{cos(a), sin(a)} * 80);
	}
	shape.push_back(double2{0, 0});
	SetShape(bodies[1], shape);

	shape.clear();
	for (int i = 0; i < 24; i++) {
		auto a = 2 * M_PI / 24 * i;
		shape.push_back(double2{cos(a), sin(a)} * 50);
	}
	shape.push_back(double2{0, 0});
	SetShape(bodies[2], shape);

	mat3 ortho;
	ortho[0][0] = 2.0f / Width;
	ortho[1][1] = 2.0f / Height;
	ortho[2][2] = 1.0f;
	ortho[2][0] = -1.0f;
	ortho[2][1] = -1.0f;

	const dvec2 gravity(0, -100); // mVm/s^2

	double time = 0;
	double energy = 0;
	for (Body& body : bodies) {
		energy += -glm::dot(body.pos, gravity) * body.mass + glm::dot(body.velocity, body.velocity) * body.mass / 2;
	}
	double energyMin = energy, energyMax = energy;

	RunEventLoop(window, [&]() {
		glClear(GL_COLOR_BUFFER_BIT);

		constexpr double dt = 0.01;
		if (gSimulate) {
			for (auto i : range(bodies.size()))
				for (auto j : range(i + 1, bodies.size()))
					Interact(bodies[i], bodies[j]);

			for (Body& body : bodies) {
				dvec2 acc = gGravity ? gravity : dvec2(0, 0);
				dvec4 s0 = dvec4(body.pos.x, body.pos.y, body.velocity.x, body.velocity.y);
				auto s1 = RungeKutta4<dvec4, double>(s0, 0, dt, [&](dvec4 s, double t) {
					return dvec4(dvec2(s.z, s.w), acc);
				});
				body.pos = dvec2(s1.x, s1.y);
				body.velocity = dvec2(s1.z, s1.w);

				constexpr double elasticity = 0.5;
				if (body.pos.y + body.box.min.y < 0) {
					double dip = -(body.pos.y + body.box.min.y);
					// remove penetration, but preserve total energy (if possible when d >= 0)
					double d = gravity.y * 2 * dip + body.velocity.y * body.velocity.y;
					body.velocity.y = (d > 0) ? elasticity * sqrt(d) : 0;
					body.pos.y = -body.box.min.y;
				}

				if (body.pos.x + body.box.min.x <= 0 && body.velocity.x < 0) {
					body.velocity.x = -body.velocity.x * elasticity;
				}
				if (body.pos.x + body.box.max.x >= Width && body.velocity.x > 0) {
					body.velocity.x = -body.velocity.x * elasticity;
				}
				if (body.pos.y + body.box.max.y >= Height && body.velocity.y > 0) {
					body.velocity.y = -body.velocity.y * elasticity;
				}

				if (gAirDrag) {
					body.velocity *= 0.999;
				}
			}
			time += dt;
		}

		double energy = 0;
		for (Body& body : bodies) {
			energy += -glm::dot(body.pos, gravity) * body.mass + glm::dot(body.velocity, body.velocity) * body.mass / 2;
		}
		minimize(energyMin, energy);
		maximize(energyMax, energy);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		for (Body& body : bodies) {
			vec2 p = body.pos;
			p.x = p.x / Width * 800;
			p.y = p.y / Height * 600;
			monaco.render(
				format("pos (%f)\nvec (%f)", body.pos, body.velocity),
				p.x, p.y, 0.3, "7FE030");
		}
		monaco.render(
			format("energy %s, spread %s, time %s", energy, energyMax - energyMin, time),
			5, 5, 0.3, "7FE030");
		glDisable(GL_BLEND);

		if (time >= 100)
			gSimulate = false;

		glUseProgram(shader);
		buffer.bind();
		for (const Body& body : bodies) {
			mat3 transform(ortho);
			transform = glm::translate(transform, vec2(body.pos));
			transform = glm::rotate(transform, float(body.angle));
			transformUniform = transform;
			buffer.write(span<const vec2>(body.shapef.data(), body.shapef.size()));
			glDrawArrays(GL_LINE_LOOP, 0, body.shapef.size());
		}
	});
	return 0;
}
