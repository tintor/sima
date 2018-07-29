#pragma once

#include <string>
#include <vector>
#include <mutex>

#define GLFW_INCLUDE_GLCOREARB
#ifndef __APPLE_CC__
	#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>

#include "glm.h"
#include "format.h"

void Error(const char* name);

int gen_buffer(GLenum target, GLsizei size, const void* data);
GLuint load_program(std::string_view name, bool geometry = false);
void load_png_texture(std::string_view filename);

class Text
{
public:
	Text();

	//void Reset(int width, int height, glm::mat4& matrix, bool down);

	template<typename ...Args>
	void Print(std::string_view fmt, Args... args) {
		PrintBuffer(format(fmt, args...));
	}
	void PrintBuffer(std::string_view s) {
		PrintAt(_tx, _ty, _ts, s);
		_ty += _tdy;
	}
	void PrintAt(float x, float y, float n, std::string_view);

	vec4 fg_color;
	vec4 bg_color;

private:
	float _tx;
	float _ty;
	float _ts;
	float _tdy;

	GLuint _positionBuffer;
	GLuint _uvBuffer;
	std::vector<GLfloat> _position_data;
	std::vector<GLfloat> _uv_data;
};

// TODO make it work for vertical monitor setup
const int ConsoleWidth = 131;
const int ConsoleHeight = 40; // TODO: magic numbers?

class Console
{
public:
	Console();
	virtual void Execute(std::string_view command) { Print("[%s]\n", command); }

	template<typename ...Args>
	void Print(std::string_view fmt, Args... args) {
		std::lock_guard lock(_mutex);
		PrintInternal(format(fmt, args...));
	}

	static bool KeyToChar(int key, int mods, char& ch);
	bool OnKey(int key, int mods);
	void Render(Text* text, float time);

	bool IsVisible() { return _visible; }
	void Show() { _visible = true; }
	void Hide() { _visible = false; }

	//json_t* save();
	//bool load(json_t* doc);

private:
	void PrintInternal(std::string_view s);

	bool _visible = false;
	std::mutex _mutex;

	char _output[ConsoleHeight][ConsoleWidth];
	int _last_line = 0;
	int _last_column = 0;

	struct Input {
		std::array<char, ConsoleWidth> buffer;
		int cursor;
	};

	int _current_input = 0;
	std::vector<Input*> _inputs;
};
