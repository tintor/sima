#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>
#include <sim/integration.h>
#include <view/font.h>
#include <view/glm.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>
#include <view/window.h>

bool gSimulate = false;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
	glfwSetWindowShouldClose(window, GL_TRUE);
	return;
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_SPACE && mods == 0) {
	gSimulate ^= 1;
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
    }
}

void scroll_callback(GLFWwindow* window, double x, double y) {}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) { glViewport(0, 0, width, height); }

struct Body {
    double radius;
    double mass;
    dvec2 pos;
    dvec2 velocity = vec2(0, 0);
    double angle = 0;
    double ang_velocity = 0;
};

int main(int argc, char** argv) {
    InitSegvHandler();

    constexpr int Width = 1200, Height = 900;
    auto window = CreateWindow({.width = Width, .height = Height, .resizeable = false});
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

    FontRenderer renderer(800, 600);
	Font timesNewRoman("Times New Roman", 48, &renderer);
	Font arial("Arial", 48, &renderer);
	Font monaco("/System/Library/Fonts/Monaco.dfont", 48, &renderer);

    std::array<vec2, 25> v;
    for (int i = 0; i < v.size() - 1; i++) {
	auto a = 2 * PI / (v.size() - 1) * i;
	v[i].x = cos(a);
	v[i].y = sin(a);
    }
    v[24] = vec2(0, 0);
    buffer.write(v);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

    vector<Body> bodies;
    bodies.resize(3);
    bodies[0] = {.mass = 2e8, .radius = 1, .pos = vec2(600, 450)};
    bodies[1] = {.mass = 600, .radius = 1, .pos = vec2(200, 450), .velocity = vec2(0, 400 / (150000000 / 30))};
    bodies[2] = {.mass = 7, .radius = 1, .pos = vec2(200 + 16 / 15, 450), .velocity = vec2(0, 16 / 15 / 400)};

    mat3 ortho;
    ortho[0][0] = 2.0f / Width;
    ortho[1][1] = 2.0f / Height;
    ortho[2][2] = 1.0f;
    ortho[2][0] = -1.0f;
    ortho[2][1] = -1.0f;

    // TODO fix
    const dvec2 gravity(0, -100);  // mVm/s^2
    double energy = 0;
    for (Body& body : bodies) {
	// TODO fix potential engery
	energy += -glm::dot(body.pos, gravity) * body.mass + glm::dot(body.velocity, body.velocity) * body.mass / 2;
    }
    double energyMin = energy, energyMax = energy;

    double time = 0;

    RunEventLoop(window, [&]() {
	glClear(GL_COLOR_BUFFER_BIT);

	constexpr double dt = 0.01;
	if (gSimulate) {
	    // F=k*M1*M2/r^2
	    // a=k*M_other/d^2
	    // TODO compute force between each pair of bodies
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
	    monaco.render(format("pos (%f)\nvec (%f)", body.pos, body.velocity), p.x, p.y, 0.3, "7FE030");
	}
	monaco.render(format("energy %s, spread %s, time %s", energy, energyMax - energyMin, time), 5, 5, 0.3,
	              "7FE030");
	glDisable(GL_BLEND);

	if (time >= 100) gSimulate = false;

	glUseProgram(shader);
	buffer.bind();
	for (const Body& body : bodies) {
	    mat3 transform(ortho);
	    transform = glm::translate(transform, vec2(body.pos));
	    transform = glm::rotate(transform, float(body.angle));
	    transform = glm::scale(transform, vec2(body.radius, body.radius));
	    transformUniform = transform;
	    glDrawArrays(GL_LINE_LOOP, 0, v.size());
	}
    });
    return 0;
}
