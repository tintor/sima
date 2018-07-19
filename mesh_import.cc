#include "mesh_import.h"
#include "util.h"
#include "file.h"

constexpr const char* FLOAT_REGEX = R"(-?\d+\.\d+(e[+-]\d+)?)";

imesh3 load_stl(std::string_view filename) {
    FileReader file(filename);
    std::string_view line;

    imesh3 mesh;
    std::cmatch m;

    if (file.readline() != "solid Default\n")
        throw error("expected [solid Default]");

    std::regex facet_regex(R"(^\s*facet\s+)");
    std::regex vertex_regex(format(R"(^\s*vertex ({}) ({}) ({})\s*$)", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (true) {
		line = file.readline();
        if (line == "endsolid Default\n")
            break;
        if (!search(line, facet_regex))
            throw error("expected [facet ...] or [endsolid Default] instead of [{}]", line);
        if (file.readline() != "    outer loop\n")
            throw error("expected [outer loop]");

        dvec3 v[3];
        for (auto i : range(3)) {
			line = file.readline();
            if (!match(line, vertex_regex, /*out*/m))
                throw error("expected [vertex ...] instead of [{}]", line);
            v[i].x = parse<double>(m[1]);
            v[i].y = parse<double>(m[3]);
            v[i].z = parse<double>(m[5]);
        }
        mesh.emplace_back(v[0], v[1], v[2]);

        if (file.readline() != "    endloop\n")
            throw error("expected [endloop]");
        if (file.readline() != "  endfacet\n")
            throw error("expected [endfacet]");
    }

    return mesh;
}

imesh3 load_ply(std::string_view filename) {
    FileReader file(filename);
    std::string_view line;

    std::vector<dvec3> vertices;
    imesh3 mesh;
    std::cmatch m;

    std::regex element_regex(R"(element (vertex|face) (\d+)\s*)");
    while (true) {
		line = file.readline();
        if (line == "")
            throw error("bad ply file header");
        if (line == "end_header\n")
            break;
        if (match(line, element_regex, /*out*/m)) {
            int s = parse<int>(m[2]);
            if (m[1].str() == "vertex")
                vertices.reserve(s);
            else
                mesh.reserve(s);
        }
    }

    std::regex vertex_regex(format(R"(^({}) ({}) ({})(\s+|$))", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (vertices.size() < vertices.capacity()) {
		line = file.readline();
        if (line == "")
            throw error("unexpected end of ply file");
        if (!search(line, vertex_regex, /*out*/m))
            throw error("expected vertex in ply file: [{}]", line);
        double a = parse<double>(m[1]);
        double b = parse<double>(m[3]);
		double c = parse<double>(m[5]);
        vertices.emplace_back(a, b, c);
    }

    std::regex face_regex(R"(3 (\d+) (\d+) (\d+)\s*)");
    while (mesh.size() < mesh.capacity()) {
        line = file.readline();
        if (line == "")
            throw error("unexpected end of ply file");
		if (!match(line, face_regex, /*out*/m))
            throw error("expected face in ply file");
        int a = parse<int>(m[1]);
        int b = parse<int>(m[3]);
		int c = parse<int>(m[5]);
        mesh.emplace_back(vertices[a], vertices[b], vertices[c]);
    }

    return mesh;
}
