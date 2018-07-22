#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <functional>
#include <stdexcept>
#include <cstdint>
#include <atomic>
#include <unistd.h>
#include <execinfo.h>

#define GLFW_INCLUDE_GLCOREARB
#ifndef __APPLE_CC__
    #include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>
#include "auto.h"
#include "rendering.h"
#include "shape.h"
#include "integration.h"
#include "timestamp.h"
#include "callstack.h"

// Dynamics and simulation
// =======================

struct Body {
    // const
    Shape shape;
    real mass;
    dmat3 inertial_tensor;
    // mutable
	dvec3 position;
	dquat orientation;
	dvec3 velocity;
	dvec3 angular_velocity;
};

class Joint {
	// ball joint (common point) 3DF
	// hinge / axel joint (common edge) 1DF
	// cylindrical joint (common line) 2DF
	// prismatic joint 1DF (two common parallel lines)
};

class BallJoint {
	dvec3 pa, pb;
};

class HingeJoint {
	dvec3 pa, da, pb, db;
	dvec3 ra, rb; // just for computing relative angle (must be unit and normal to da / db)
};

// Open chain articulated only
// TODO describe children with relative coordinates, but allow computation of absolute ones
class ArticulatedBody {
	Body m_body;
	std::vector<std::pair<ArticulatedBody, Joint>> m_children;
};

// TODO update to use custom integrator
void update(Body& body, double dt) {
    dvec3 bias_velocity, bias_angular_velocity;
    body.position += dt * (body.velocity + bias_velocity);
    // dq/dt = w*q/2  =>  q' = q + (w*q)*(dt/2)
    //body.angular_velocity * body.orientation;
    //body.orientation = glm::normalize(body.orientation + (body.angular_velocity + bias_angular_velocity) * body.orientation * (dt / 2));

    dvec3 torque, force;
    body.velocity += force * (dt / body.mass);
    auto invI = glm::inverse(body.inertial_tensor);
    auto I = body.inertial_tensor;
    body.angular_velocity += invI * (torque - glm::cross(body.angular_velocity, I * body.angular_velocity)) * dt;
}


// =====================

Text* text = nullptr;

struct FpvCamera {
    vec3 position;
    quat orientation;
};
FpvCamera camera;

glm::vec3 g_position(0,0,0);
float g_yaw=0, g_pitch=0;
glm::mat4 g_orientation;

void on_key(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
		glfwSetWindowShouldClose(window, GL_TRUE);
		return;
	}
	// TODO:
	// W and S, forward / backward
	// Q and E, roll
	// update g_orientation of mouse move
}

void on_mouse_button(GLFWwindow* window, int button, int action, int mods) {
	if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
	}
}

void on_scroll(GLFWwindow* window, double x, double y) {
}

// =====================
// OpenGL

int render_width = 0, render_height = 0;

void model_init(GLFWwindow* window) {
    glfwSetKeyCallback(window, on_key);
    glfwSetMouseButtonCallback(window, on_mouse_button);
    glfwSetScrollCallback(window, on_scroll);
}

glm::mat4 perspective, perspective_rotation;

void render_init() {
	fprintf(stderr, "OpenGL version: [%s]\n", glGetString(GL_VERSION));
	glEnable(GL_CULL_FACE);
	glClearColor(0.0, 0.5, 0.0, 1.0);
	glViewport(0, 0, render_width, render_height);

	::perspective = glm::perspective<float>(M_PI / 180 * 90, render_width / (float)render_height, 0.03, 1000);
	perspective_rotation = glm::rotate<float>(::perspective, 0, glm::vec3(1, 0, 0));
	perspective_rotation = glm::rotate<float>(perspective_rotation, 0, glm::vec3(0, 1, 0));
	perspective_rotation = glm::rotate<float>(perspective_rotation, float(M_PI / 2), glm::vec3(-1, 0, 0));

	text = new Text;
	text->fg_color = vec4(1, 1, 1, 1);
	text->bg_color = vec4(0, 0, 0, 1);
}

struct Triangle {
    vec3 vertex[3];
    vec3 color;
};

struct Line {
    vec3 vectex[2];
    vec3 color;
};

void render_line(vec3 vertex_a, vec3 vertex_b, vec3 color) {

}

void render_world(const glm::mat4& matrix) {
    glClear(GL_COLOR_BUFFER_BIT /*| GL_DEPTH_BUFFER_BIT*/);
    //glEnable(GL_DEPTH_TEST);
    // TODO floor checkerbox
    // TODO two boxes in space
    //glDisable(GL_DEPTH_TEST);
    render_line(vec3(0,0,0), vec3(3,0,0), vec3(1,1,1));
    render_line(vec3(0,0,0), vec3(0,2,0), vec3(1,1,1));
    render_line(vec3(0,0,0), vec3(0,0,1), vec3(1,1,1));
}

void render_gui() {
	glm::mat4 matrix = glm::ortho<float>(0, render_width, 0, render_height, -1, 1);
	text->Reset(render_width, render_height, matrix, true);
	text->Print("Hello world!");
}

bool last_cursor_init = false;
double last_cursor_x, last_cursor_y;

void turn(float dx, float dy) {
        g_yaw += dx;
        g_pitch += dy;
        if (g_pitch > M_PI / 2 * 0.999) g_pitch = M_PI / 2 * 0.999;
        if (g_pitch < -M_PI / 2 * 0.999) g_pitch = -M_PI / 2 * 0.999;
        g_orientation = glm::rotate(glm::rotate(glm::mat4(), -g_yaw, glm::vec3(0, 0, 1)), -g_pitch, glm::vec3(1, 0, 0));
}

void model_orientation(GLFWwindow* window) {
        double cursor_x, cursor_y;
        glfwGetCursorPos(window, &cursor_x, &cursor_y);
        if (!last_cursor_init) {
                last_cursor_init = true;
                last_cursor_x = cursor_x;
                last_cursor_y = cursor_y;
        }
        if (cursor_x != last_cursor_x || cursor_y != last_cursor_y) {
                turn((cursor_x - last_cursor_x) / 150, (cursor_y - last_cursor_y) / 150);
                last_cursor_x = cursor_x;
                last_cursor_y = cursor_y;
                perspective_rotation = glm::rotate(::perspective, g_pitch, glm::vec3(1, 0, 0));
                perspective_rotation = glm::rotate(perspective_rotation, g_yaw, glm::vec3(0, 1, 0));
                perspective_rotation = glm::rotate(perspective_rotation, float(M_PI / 2), glm::vec3(-1, 0, 0));
        }
}

void OnError(int error, const char* message) {
	fprintf(stderr, "GLFW error %d: %s\n", error, message);
}

GLFWwindow* create_window() {
	glfwSetErrorCallback(OnError);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
	GLFWwindow* window = glfwCreateWindow(mode->width * 2, mode->height * 2, "Sima", glfwGetPrimaryMonitor(), NULL);
	if (!window)
		return nullptr;
	glfwMakeContextCurrent(window);
	glfwSwapInterval(0/*VSYNC*/);
	//glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwGetFramebufferSize(window, &render_width, &render_height);
	return window;
}

void sigsegv_handler(int sig) {
	fprintf(stderr, "Error: signal %d:\n", sig);
	Callstack stack;
	std::string s;
	stack.write(s, {"sigsegv_handler(int)", "_sigtramp"});
	fputs(s.c_str(), stderr);
	exit(1);
}

// solid angle between triangle and origin
real solid_angle(dvec3 A, dvec3 B, dvec3 C) {
    real y = glm::dot(A, glm::cross(B, C));
    real a = sqrt(squared(A)), b = sqrt(squared(B)), c = sqrt(squared(C));
    real x = a * b * c + c * glm::dot(A, B) + b * glm::dot(A, C) + a * glm::dot(B, C);
    return 2 * std::atan2(y, x);
}

int main(int argc, char** argv) {
	void sigsegv_handler(int sig);
	signal(SIGSEGV, sigsegv_handler);

    Timestamp::init();

    if (!glfwInit())
        return -1;

    /*try {
        Mesh3d mm = load_stl("models/bunny.stl");
        std::vector<dvec3> vertices;
        FOR_EACH(f, mm)
            FOR(i, 3)
                vertices.push_back(f[i]);
        Mesh3d ch = build_convex_hull(vertices);
        std::cout << "IsValid " << static_cast<int>(is_valid(mm)) << std::endl;
        std::cout << "Volume " << volume(mm) << std::endl;
        std::cout << "CenterOfMass " << center_of_mass(mm) << std::endl;
        std::cout << "IsConvex " << is_convex(mm) << std::endl;

        std::cout << "IsValid " << static_cast<int>(is_valid(ch)) << std::endl;
        std::cout << "Volume " << volume(ch) << std::endl;
        std::cout << "CenterOfMass " << center_of_mass(ch) << std::endl;
        std::cout << "IsConvex " << is_convex(ch) << std::endl;

        std::default_random_engine rnd(0);

        Timestamp ta;
        SolidBSPTree tree(mm, 100000, rnd);
        std::cout << ta.elapsed_ms() << std::endl;
    } catch (std::runtime_error& e) {
        std::cout << "std::runtime_error " << e.what() << std::endl;
    }
    return 0;*/

	GLFWwindow* window = create_window();
	if (!window)
		return -1;
	model_init(window);
	render_init();

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		//model_frame(window, frame_ms);
		glm::mat4 matrix = glm::translate(perspective_rotation, g_position);
		//Frustum frustum(matrix);
		//g_player.cpos = glm::ivec3(glm::floor(g_player.position)) >> ChunkSizeBits;*/

        render_world(matrix);
		render_gui();

		glfwSwapBuffers(window);
	}

	glfwTerminate();
	_exit(0); // exit(0) is not enough
	return 0;
}
