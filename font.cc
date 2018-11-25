#include "font.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <ft2build.h>
#include FT_FREETYPE_H

#include "shader.h"
#include "auto.h"
#include "format.h"

struct Renderer {
	uint VAO, VBO;
	Shader shader;
	int textColorLocation;

	Renderer();
};

Renderer::Renderer() : shader(R"END(
		#version 330 core
		layout (location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>
		out vec2 TexCoords;

		uniform mat4 projection;

		void main() {
		    gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);
		    TexCoords = vertex.zw;
		}

		#version 330 core
		in vec2 TexCoords;
		out vec4 color;

		uniform sampler2D text;
		uniform vec3 textColor;

		void main() {
		    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);
		    color = vec4(textColor, 1.0) * sampled;
		}
	)END") {
	glm::mat4 projection = glm::ortho(0.0, 800.0, 0.0, 600.0);
	shader.use();
	glUniformMatrix4fv(glGetUniformLocation(shader, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
	textColorLocation = glGetUniformLocation(shader, "textColor");

	// Disable byte-alignment restriction
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// Configure VAO/VBO for texture quads
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, nullptr, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

static Renderer* s_renderer = nullptr;

Font::Font(string_view name, int resolution) {
	if (s_renderer == nullptr) {
		s_renderer = new Renderer();
	}

    FT_Library ft;
    if (FT_Init_FreeType(&ft)) {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
		exit(0);
	}
	ON_SCOPE_EXIT(FT_Done_FreeType(ft));

    // Load font as face
	// TODO do search in both: /System/Library/Fonts and /Library/Fonts with various extensions
    const char* fmt = name.find('.') != string_view::npos ? "%s" : "/Library/Fonts/%s.ttf";
    FT_Face face;
    if (FT_New_Face(ft, format(fmt, name).c_str(), 0, &face)) {
        std::cout << "ERROR::FREETYPE: Failed to load font: " << format(fmt, name) << std::endl;
		exit(0);
	}
	ON_SCOPE_EXIT(FT_Done_Face(face));

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, resolution);

	ON_SCOPE_EXIT(glBindTexture(GL_TEXTURE_2D, 0));
    for (int c = 0; c < m_characters.size(); c++) {
        // Load character glyph
        if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
			printf("ERROR::FREETYTPE: Failed to load Glyph '%c' (%d)\n", c, c);
			exit(0);
        }

		GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer
        );

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		m_characters[c] = {
            texture,
            face->glyph->bitmap.width, face->glyph->bitmap.rows,
            face->glyph->bitmap_left, face->glyph->bitmap_top,
            uint(face->glyph->advance.x)
        };
    }

	m_max_size_y = 0;
	for (const Character& ch : m_characters) {
		m_max_size_y = std::max(m_max_size_y, ch.size_y);
	}
}

Font::~Font() {
	// TODO delete textures
}

void Font::render(string_view text, double x, double y, double scale, Color color) {
	const auto orig_x = x;
    s_renderer->shader.use();
    glUniform3f(s_renderer->textColorLocation, color.r, color.g, color.b);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(s_renderer->VAO);

    for (char c : text) {
		if (c == '\n') {
			x = orig_x;
			y -= m_max_size_y * scale * 1.2;
			continue;
		}

        Character ch = m_characters[c];
        GLfloat xpos = x + ch.bearing_x * scale;
        GLfloat ypos = y - (int(ch.size_y) - int(ch.bearing_y)) * scale;

        GLfloat w = ch.size_x * scale;
        GLfloat h = ch.size_y * scale;
        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }
        };
        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.texture);
        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, s_renderer->VBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)
		// Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
        x += (ch.advance >> 6) * scale;
    }
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
}
