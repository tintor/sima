#pragma once

#include <glad/glad.h>
#include "std.h"
#include "auto.h"

class Shader
{
public:
    bool load(string_view vertexCode, string_view fragmentCode) {
        const char* vShaderCode = vertexCode.data();
        const char* fShaderCode = fragmentCode.data();
		int vertexLength = vertexCode.size();
		int fragmentLength = fragmentCode.size();

        // vertex shader
        uint vertex = glCreateShader(GL_VERTEX_SHADER);
        ON_SCOPE_EXIT(glDeleteShader(vertex));
        glShaderSource(vertex, 1, &vShaderCode, &vertexLength);
        glCompileShader(vertex);
        if (!checkCompileErrors(vertex, "VERTEX"))
			return false;
        // fragment Shader
        uint fragment = glCreateShader(GL_FRAGMENT_SHADER);
        ON_SCOPE_EXIT(glDeleteShader(fragment));
        glShaderSource(fragment, 1, &fShaderCode, &fragmentLength);
        glCompileShader(fragment);
        if (!checkCompileErrors(fragment, "FRAGMENT"))
			return false;
        // shader Program
        m_id = glCreateProgram();
        glAttachShader(m_id, vertex);
        glAttachShader(m_id, fragment);
        glLinkProgram(m_id);
        checkCompileErrors(m_id, "PROGRAM");
		return true;
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
    bool checkCompileErrors(uint shader, string_view type) {
        int success;
        char infoLog[1024];
        if (type != "PROGRAM") {
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
            if (!success) {
                glGetShaderInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n" << infoLog << std::endl;
				return false;
            }
        } else {
            glGetProgramiv(shader, GL_LINK_STATUS, &success);
            if (!success) {
                glGetProgramInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n" << infoLog << std::endl;
				return false;
            }
        }
		return true;
    }

	uint m_id = 0;
};
