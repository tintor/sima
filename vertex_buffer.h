#pragma once

#include "std.h"
#include "glm.h"

#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>

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
