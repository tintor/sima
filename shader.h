#pragma once

#include "auto.h"
#include "format.h"

#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>

class Shader
{
public:
	Shader(string_view source) {
		size_t n = source.rfind("#version ");
		string_view vert = source.substr(0, n);
		string_view frag = source.substr(n);

        uint vertex = compile(GL_VERTEX_SHADER, vert);
        uint fragment = compile(GL_FRAGMENT_SHADER, frag);
        ON_SCOPE_EXIT(glDeleteShader(vertex));
        ON_SCOPE_EXIT(glDeleteShader(fragment));

        m_id = glCreateProgram();
        glAttachShader(m_id, vertex);
        glAttachShader(m_id, fragment);
        glLinkProgram(m_id);

        int success;
        glGetProgramiv(m_id, GL_LINK_STATUS, &success);
        if (!success) {
        	array<char, 1024> infoLog;
            glGetProgramInfoLog(m_id, infoLog.size(), nullptr, infoLog.data());
			print("PROGRAM LINKING ERROR:\n");
			print("%s\n", infoLog);
			exit(0);
        }
	}

	operator uint() const { return m_id; }

	void use() {
        glUseProgram(m_id);
    }

	void setBool(const char* name, bool value) const {
        glUniform1i(glGetUniformLocation(m_id, name), (int)value);
    }

    void setInt(const char* name, int value) const {
        glUniform1i(glGetUniformLocation(m_id, name), value);
    }

    void setFloat(const char* name, float value) const {
        glUniform1f(glGetUniformLocation(m_id, name), value);
    }

private:
	uint compile(int type, string_view source) {
        const char* ptr = source.data();
		int length = source.size();

        uint shader = glCreateShader(type);
        glShaderSource(shader, 1, &ptr, &length);
        glCompileShader(shader);
        int success;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
        	array<char, 1024> infoLog;
            glGetShaderInfoLog(shader, infoLog.size(), nullptr, infoLog.data());
			print("%s SHADER COMPILATION ERROR:\n", (type == GL_VERTEX_SHADER) ? "VERTEX" : "FRAGMENT");
			print("%s\n", infoLog);
			exit(0);
        }
		return shader;
	}

private:
	uint m_id = 0;
};
