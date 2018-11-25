#pragma once
#include "std.h"
#include "color.h"

// TODO: ANTI-ALIAS: Valve's signed distance fields:
//       - https://steamcdn-a.akamaihd.net/apps/valve/2007/SIGGRAPH2007_AlphaTestedMagnification.pdf
// TODO: custom window sizes
// TODO: unicode characters
// TODO: x and y positioning: left, center, right
// TODO: 3d text positioning
// TODO: text rotation
// TODO: compute text size

// TODO: markup language to change text style and color
// TODO: text console

class Font {
public:
	Font(string_view name, int resolution = 48);
	~Font();
	// range for (x, y) is always (800, 600)
	void render(string_view text, double x, double y, double scale, Color color);

private:
	int m_max_size_y = 0;

	uint m_texture = 0;

	struct Character {
		float u0;
		float u1;
		float v0;
		float v1;

		int size_x, size_y; // Size of glyph
		int bearing_x, bearing_y; // Offset from baseline to left/top of glyph
	    int advance;    // Horizontal offset to advance to next glyph
	};

	array<Character, 128> m_characters;
};
