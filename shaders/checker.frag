#version 150 core

in vec2 fragment_uv;
out vec4 color;

void main() {
    int f = int(fragment_uv.x) ^ int(fragment_uv.y);
    color = (f & 1) == 0 ? vec4(0, 1, 0, 0) : vec4(0, 0.5, 0, 0);
}
