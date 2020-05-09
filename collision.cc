#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>
#include <geom/classify.h>
#include <geom/pose.h>
#include <geom/primitives.h>
#include <geom/properties.h>
#include <geom/triangle.h>
#include <view/font.h>
#include <view/glm.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>
#include <view/window.h>

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    const char* key_name = glfwGetKeyName(key, 0);
    print("key_callback [%s] key:%s scancode:%s action:%s mods:%s\n", key_name, key, scancode, action, mods);

    if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
	glfwSetWindowShouldClose(window, GL_TRUE);
	return;
    }
}

vector<polygon2> shapes;
double2 prev_cursor;
double2 cursor;
bool dragging = false;
int drag_shape = -1;
int drag_index = -1;
double2 drag_start;

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (!dragging && action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
	double min_d = 1e100;
	for (size_t i = 0; i < shapes.size(); i++) {
	    for (size_t j = 0; j < shapes[i].size(); j++) {
		double d = length(shapes[i][j] - cursor);
		if (d < min_d && d <= 5) {
		    min_d = d;
		    drag_shape = i;
		    drag_index = j;
		    dragging = true;
		    prev_cursor = cursor;
		    drag_start = cursor;
		}
	    }
	}
    }

    if (dragging && action == GLFW_RELEASE && button == GLFW_MOUSE_BUTTON_LEFT) {
	dragging = false;
    }
}

void scroll_callback(GLFWwindow* window, double x, double y) {}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) { glViewport(0, 0, width, height); }

constexpr int Width = 1600, Height = 1200;

double distance(const polygon2& a, const polygon2& b) {
    double d = 1e100;
    for (const double2& p : a)
	for (auto e : Edges(b)) d = std::min(distance(p, e), d);
    for (const double2& p : b)
	for (auto e : Edges(a)) d = std::min(distance(p, e), d);
    return d;
}

int main(int argc, char** argv) {
    InitSegvHandler();

    auto window = CreateWindow({.width = Width, .height = Height, .resizeable = false});
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);

    glClearColor(0.0, 0.5, 0.0, 1.0);
    glPointSize(10);

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

    VertexBuffer_double2 buffer(25);

    Font timesNewRoman("Times New Roman");
    Font arial("Arial");
    Font monaco("/System/Library/Fonts/Monaco.dfont");

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

    // Shape 1
    polygon2 s;
    s.push_back(double2{50, 250});
    s.push_back(double2{250, 250});
    s.push_back(double2{250, 150});
    s.push_back(double2{150, 150});
    s.push_back(double2{150, 50});
    s.push_back(double2{50, 50});
    shapes.push_back(s);

    // Shape 2
    s.clear();
    for (int i = 0; i < 6; i++) {
	auto a = 2 * PI / 6 * i;
	s.push_back(double2{cos(a), sin(a)} * 100 + double2{500, 500});
    }
    shapes.push_back(s);

    mat3 ortho;
    ortho[0][0] = 2.0f / Width;
    ortho[1][1] = 2.0f / Height;
    ortho[2][2] = 1.0f;
    ortho[2][0] = -1.0f;
    ortho[2][1] = -1.0f;

    RunEventLoop(window, [&]() {
	prev_cursor = cursor;
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	cursor = double2{x, Height - y};
	if (dragging) {
	    double2 a = drag_start;
	    double2 b = cursor;

	    shapes[drag_shape][drag_index] = b;
	    int c = Classify(shapes[0], shapes[1]);
	    if (c == -1) {
		double aa = 0;
		double bb = 1e100;

		double2 valid = a;
		double vv = 0;

		for (int i = 0; i < 40; i++) {
		    double2 mid = (a + b) / 2;
		    double mm = (aa + bb) / 2;
		    shapes[drag_shape][drag_index] = mid;
		    int c = Classify(shapes[0], shapes[1]);
		    if (c != -1 && mm > vv) {
			vv = mm;
			valid = mid;
		    }
		    // print("%d %.12f [%d]\n", i, distance(shapes[0], shapes[1]), c);
		    if (c == -1) {
			b = mid;
			bb = mm;
		    } else {
			a = mid;
			aa = mm;
		    }
		}
		b = valid;
	    }
	    shapes[drag_shape][drag_index] = b;
	}

	glClear(GL_COLOR_BUFFER_BIT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	for (size_t a = 0; a < shapes.size(); a++) {
	    glEnable(GL_BLEND);
	    for (size_t b = a + 1; b < shapes.size(); b++) {
		vector<IContact2> contacts;
		int c = Classify(shapes[a], shapes[b], &contacts);
		string s;
		for (const auto& con : contacts) {
		    s += format("\nnormal:(%.6f %.6f) sa:(%.6f %.6f)", con.normal.x, con.normal.y, con.sa.x, con.sa.y);
		}
		monaco.render(format("shape %d %d -> %d | dist %.12f%s", a, b, c, distance(shapes[a], shapes[b]), s),
		              100, 100, 0.2, "FFFFFF");
	    }
	    glDisable(GL_BLEND);

	    glUseProgram(shader);
	    buffer.bind();
	    mat3 transform(ortho);
	    transformUniform = transform;
	    buffer.write(shapes[a]);
	    glDrawArrays(GL_LINE_LOOP, 0, shapes[a].size());
	    glDrawArrays(GL_POINTS, 0, shapes[a].size());
	}
    });
    return 0;
}
