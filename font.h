#pragma once
#include "std.h"
#include "color.h"

// TODO: ANTI-ALIAS: Valve's signed distance fields
// TODO: PERF: use texture-per-font instead of texture-per-character
// TODO: custom window sizes
// TODO: unicode characters
// TODO: x and y positioning: left, center, right
// TODO: 3d text positioning
// TODO: use hex format for color
// TODO: text rotation
// TODO: compute text size

class Font {
public:
	Font(string_view name, int resolution = 48);
	~Font();
	// range for (x, y) is always (800, 600)
	void render(string_view text, double x, double y, double scale, Color color);

private:
	uint m_max_size_y = 0;
	struct Character {
	   	uint texture;   // ID handle of the glyph texture
		uint size_x, size_y; // Size of glyph
		int bearing_x, bearing_y; // Offset from baseline to left/top of glyph
	    uint advance;    // Horizontal offset to advance to next glyph
	};

	array<Character, 128> m_characters;
};
