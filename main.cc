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

#include "solid_bsp_tree.h"

// Dynamics and simulation
// =======================

struct Body {
    // const
    Shape shape;
    double mass;
    //dmat3 inertial_tensor;
    // mutable
	double3 position;
	double4 orientation;
	double3 velocity;
	double3 angular_velocity;
};

class Joint {
	// ball joint (common point) 3DF
	// hinge / axel joint (common edge) 1DF
	// cylindrical joint (common line) 2DF
	// prismatic joint 1DF (two common parallel lines)
};

class BallJoint {
	double3 pa, pb;
};

class HingeJoint {
	double3 pa, da, pb, db;
	double3 ra, rb; // just for computing relative angle (must be unit and normal to da / db)
};

// Open chain articulated only
// TODO describe children with relative coordinates, but allow computation of absolute ones
class ArticulatedBody {
	Body m_body;
	vector<pair<ArticulatedBody, Joint>> m_children;
};

// TODO update to use custom integrator
void update(Body& body, double dt) {
    double3 bias_velocity, bias_angular_velocity;
    body.position += dt * (body.velocity + bias_velocity);
    // dq/dt = w*q/2  =>  q' = q + (w*q)*(dt/2)
    //body.angular_velocity * body.orientation;
    //body.orientation = glm::normalize(body.orientation + (body.angular_velocity + bias_angular_velocity) * body.orientation * (dt / 2));

    double3 torque, force;
    body.velocity += force * (dt / body.mass);
    //auto invI = glm::inverse(body.inertial_tensor);
    //auto I = body.inertial_tensor;
    // body.angular_velocity += invI * (torque - glm::cross(body.angular_velocity, I * body.angular_velocity)) * dt;
}


// =====================

Text* text = nullptr;

struct FpvCamera {
    double3 position;
    double4 orientation;
};
FpvCamera camera;

double3 g_position{0,0,0};
float g_yaw=0, g_pitch=0;
//glm::mat4 g_orientation;

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

//glm::mat4 perspective, perspective_rotation;

/*mat44 perspective(float fovy, float aspect, float zNear, float zFar) {
	T const tanHalfFovy = tan(fovy / static_cast<T>(2));
	mat44 Result(static_cast<T>(0));
	Result[0][0] = static_cast<T>(1) / (aspect * tanHalfFovy);
	Result[1][1] = static_cast<T>(1) / (tanHalfFovy);
	Result[2][2] = - (zFar + zNear) / (zFar - zNear);
	Result[2][3] = - static_cast<T>(1);
	Result[3][2] = - (static_cast<T>(2) * zFar * zNear) / (zFar - zNear);
	return Result;
}*/

/*mat44 rotate(mat44 const& m, float angle, double3 const& v) {
	T const a = angle;
	T const c = cos(a);
	T const s = sin(a);

	vec<3, T, Q> axis(normalize(v));
	vec<3, T, Q> temp((T(1) - c) * axis);

	mat<4, 4, T, Q> Rotate;
	Rotate[0][0] = c + temp[0] * axis[0];
	Rotate[0][1] = temp[0] * axis[1] + s * axis[2];
	Rotate[0][2] = temp[0] * axis[2] - s * axis[1];

	Rotate[1][0] = temp[1] * axis[0] - s * axis[2];
	Rotate[1][1] = c + temp[1] * axis[1];
	Rotate[1][2] = temp[1] * axis[2] + s * axis[0];

	Rotate[2][0] = temp[2] * axis[0] + s * axis[1];
	Rotate[2][1] = temp[2] * axis[1] - s * axis[0];
	Rotate[2][2] = c + temp[2] * axis[2];

	mat44 Result;
	Result[0] = m[0] * Rotate[0][0] + m[1] * Rotate[0][1] + m[2] * Rotate[0][2];
	Result[1] = m[0] * Rotate[1][0] + m[1] * Rotate[1][1] + m[2] * Rotate[1][2];
	Result[2] = m[0] * Rotate[2][0] + m[1] * Rotate[2][1] + m[2] * Rotate[2][2];
	Result[3] = m[3];
	return Result;
}
*/

void render_init() {
	fprintf(stderr, "OpenGL version: [%s]\n", glGetString(GL_VERSION));
	glEnable(GL_CULL_FACE);
	glClearColor(0.0, 0.5, 0.0, 1.0);
	glViewport(0, 0, render_width, render_height);

	/*::perspective = glm::perspective<float>(M_PI / 180 * 90, render_width / (float)render_height, 0.03, 1000);
	perspective_rotation = glm::rotate<float>(::perspective, 0, glm::double3(1, 0, 0));
	perspective_rotation = glm::rotate<float>(perspective_rotation, 0, glm::double3(0, 1, 0));
	perspective_rotation = glm::rotate<float>(perspective_rotation, float(M_PI / 2), glm::double3(-1, 0, 0));*/

	text = new Text;
	text->fg_color = float4{1, 1, 1, 1};
	text->bg_color = float4{0, 0, 0, 1};
}

struct Triangle {
    double3 vertex[3];
    double3 color;
};

struct Line {
    double3 vectex[2];
    double3 color;
};

void render_line(double3 vertex_a, double3 vertex_b, double3 color) {

}

#ifdef xxx
void render_world(const glm::mat4& matrix) {
    glClear(GL_COLOR_BUFFER_BIT /*| GL_DEPTH_BUFFER_BIT*/);
    //glEnable(GL_DEPTH_TEST);
    // TODO floor checkerbox
    // TODO two boxes in space
    //glDisable(GL_DEPTH_TEST);
    render_line(double3{0,0,0}, double3{3,0,0}, double3{1,1,1});
    render_line(double3{0,0,0}, double3{0,2,0}, double3{1,1,1});
    render_line(double3{0,0,0}, double3{0,0,1}, double3{1,1,1});
}
#endif

void render_gui() {
	/*glm::mat4 matrix = glm::ortho<float>(0, render_width, 0, render_height, -1, 1);
	text->Reset(render_width, render_height, matrix, true);
	text->Print("Hello world!");*/
}

bool last_cursor_init = false;
double last_cursor_x, last_cursor_y;

void turn(float dx, float dy) {
        g_yaw += dx;
        g_pitch += dy;
        if (g_pitch > M_PI / 2 * 0.999) g_pitch = M_PI / 2 * 0.999;
        if (g_pitch < -M_PI / 2 * 0.999) g_pitch = -M_PI / 2 * 0.999;
        //g_orientation = glm::rotate(glm::rotate(glm::mat4(), -g_yaw, glm::double3(0, 0, 1)), -g_pitch, glm::double3(1, 0, 0));
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
                //perspective_rotation = glm::rotate(::perspective, g_pitch, glm::double3(1, 0, 0));
                //perspective_rotation = glm::rotate(perspective_rotation, g_yaw, glm::double3(0, 1, 0));
                //perspective_rotation = glm::rotate(perspective_rotation, float(M_PI / 2), glm::double3(-1, 0, 0));
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

int main(int argc, char** argv) {
	void sigsegv_handler(int sig);
	signal(SIGSEGV, sigsegv_handler);

    Timestamp::init();

    if (!glfwInit())
        return -1;

	GLFWwindow* window = create_window();
	if (!window)
		return -1;
	model_init(window);
	render_init();

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		//model_frame(window, frame_ms);
		//glm::mat4 matrix = glm::translate(perspective_rotation, g_position);
		//Frustum frustum(matrix);
		//g_player.cpos = glm::idouble3(glm::floor(g_player.position)) >> ChunkSizeBits;*/

        //render_world(matrix);
		render_gui();

		glfwSwapBuffers(window);
	}

	glfwTerminate();
	_exit(0); // exit(0) is not enough
	return 0;
}
