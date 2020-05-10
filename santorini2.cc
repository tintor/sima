#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>
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

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
    }
}

void scroll_callback(GLFWwindow* window, double x, double y) {}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) { glViewport(0, 0, width, height); }

constexpr int Width = 1000, Height = 1000;

enum class Figure : char { None, Dome, Player1, Player2 };

struct Cell {
    char level = 0;
    Figure figure = Figure::None;
};

struct Board {
    std::array<std::array<Cell, 5>, 5> cell;
};

struct View : public VertexBuffer_vec2_rgba {
    Shader shader;
    Uniform_mat4 transform;

    FontRenderer font_renderer;
    Font mono;
    glm::mat4 ortho;

    View()
        : VertexBuffer_vec2_rgba(25),
          shader(R"END(
    		#version 330 core
    		layout (location = 0) in vec2 pos;
    		layout (location = 1) in vec4 rgba;
    		out vec4 pixel_rgba;
    		uniform mat4 transform;

            void main() {
    		    vec4 p = transform * vec4(pos, 0.0, 1.0);
		        gl_Position = p;
                pixel_rgba = rgba;
	    	}

		    #version 330 core
    		in vec4 pixel_rgba;
    		out vec4 color;

    		void main() {
    		    color = pixel_rgba;
    		}
        	)END"),
          font_renderer(1000, 1000),
          mono("JetBrainsMono-Medium.ttf", 48, &font_renderer),
          transform("transform") {}
};

void Render(const Board& board, View& view) {
    glUseProgram(view.shader);
    view.bind();
    view.transform = view.ortho;

    // Frame
    for (int i = 0; i <= 5; i++) {
	const uint64_t color = 0xFF808080;
	view.add({100, 100 + i * 160}, color);
	view.add({100 + 800, 100 + i * 160}, color);

	view.add({100 + i * 160, 100}, color);
	view.add({100 + i * 160, 100 + 800}, color);
    }
    view.draw(GL_LINES);

    for (int y = 0; y < 5; y++) {
	for (int x = 0; x < 5; x++) {
	    const double s = 160;
	    const double px = 100 + x * s + s / 2;
	    const double py = 100 + y * s + s / 2;
	    const auto cell = board.cell[x][y];

	    // Towers
	    if (cell.level >= 1) {
		const uint64_t color = 0xFFFFFFFF;
		double a = 70;
		view.add({px - a, py - a}, color);
		view.add({px + a, py - a}, color);
		view.add({px + a, py + a}, color);
		view.add({px - a, py + a}, color);
		view.draw(GL_LINE_LOOP);
	    }
	    if (cell.level >= 2) {
		const uint64_t color = 0xFFFFFFFF;
		double a = 60;
		view.add({px - a, py - a}, color);
		view.add({px + a, py - a}, color);
		view.add({px + a, py + a}, color);
		view.add({px - a, py + a}, color);
		view.draw(GL_LINE_LOOP);
	    }
	    if (cell.level >= 3) {
		const uint64_t color = 0xFFFFFFFF;
		double a = 50;
		view.add({px - a, py - a}, color);
		view.add({px + a, py - a}, color);
		view.add({px + a, py + a}, color);
		view.add({px - a, py + a}, color);
		view.draw(GL_LINE_LOOP);
	    }

	    // Domes
	    if (cell.figure == Figure::Dome) {
		const uint64_t color = 0xFFFF0000;
		double a = 40;
		for (int i = 0; i < 16; i++) {
		    double k = i * 2 * M_PI / 16;
		    view.add({px + cos(k) * a, py + sin(k) * a}, color);
		}
		view.draw(GL_TRIANGLE_FAN);
	    }

	    // Player pieces
	    if (cell.figure == Figure::Player1 || cell.figure == Figure::Player2) {
		const uint64_t color = (cell.figure == Figure::Player1) ? 0xFF00FFFF : 0xFF0000FF;
		double a = 40;
		view.add({px - a, py - a}, color);
		view.add({px + a, py - a}, color);
		view.add({px + a, py + a}, color);
		view.add({px - a, py + a}, color);
		/*for (int i = 0; i < 16; i++) {
		    double k = i * 2 * M_PI / 16;
		    view.add({px + cos(k) * a, py + sin(k) * a}, color);
		}*/
		view.draw(GL_TRIANGLE_FAN);
	    }
	}
    }
}

int main(int argc, char** argv) {
    InitSegvHandler();

    auto window = CreateWindow({.width = Width, .height = Height, .resizeable = false});
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glViewport(0, 0, Width, Height);

    glClearColor(0.0, 0.0, 0.0, 1.0);

    View view;
    view.ortho = glm::ortho(0.0, double(Width), 0.0, double(Height));

    Board board;
    board.cell[2][1].level = 3;
    board.cell[2][1].figure = Figure::Dome;
    board.cell[0][4].figure = Figure::Player1;
    board.cell[1][4].figure = Figure::Player2;
    board.cell[0][3].figure = Figure::Dome;

    RunEventLoop(window, [&]() {
	glClear(GL_COLOR_BUFFER_BIT);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	view.mono.render(format("Hello world,\nMarko!"), 500, 500, 13.0 / 48, "FFFFFF");
	glDisable(GL_BLEND);

	Render(board, view);
    });
    return 0;
}
