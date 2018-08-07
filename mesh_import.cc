#include "mesh_import.h"
#include "util.h"
#include "file.h"
#include "range.h"

constexpr const char* FLOAT_REGEX = R"(-?\d+\.\d+(e[+-]\d+)?)";

imesh3 load_stl(std::string_view filename, double scale) {
    FileReader file(filename);
    std::string_view line;

    imesh3 mesh;
    std::cmatch m;

    if (file.readline() != "solid Default\n")
        THROW(runtime_error, "expected [solid Default]");

    std::regex facet_regex(R"(^\s*facet\s+)");
    std::regex vertex_regex(format(R"(^\s*vertex (%s) (%s) (%s)\s*$)", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (true) {
		line = file.readline();
        if (line == "endsolid Default\n")
            break;
        if (!search(line, facet_regex))
            THROW(runtime_error, "expected [facet ...] or [endsolid Default] instead of [%s]", line);
        if (file.readline() != "    outer loop\n")
            THROW(runtime_error, "expected [outer loop]");

        dvec3 v[3];
        for (auto i : range(3)) {
			line = file.readline();
            if (!match(line, vertex_regex, /*out*/m))
                THROW(runtime_error, "expected [vertex ...] instead of [%s]", line);
            v[i].x = parse<double>(m[1]);
            v[i].y = parse<double>(m[3]);
            v[i].z = parse<double>(m[5]);
			v[i] *= scale;
        }
        mesh.emplace_back(vconvert(v[0], ivec3), vconvert(v[1], ivec3), vconvert(v[2], ivec3));

        if (file.readline() != "    endloop\n")
            THROW(runtime_error, "expected [endloop]");
        if (file.readline() != "  endfacet\n")
            THROW(runtime_error, "expected [endfacet]");
    }

    return mesh;
}

imesh3 load_ply(std::string_view filename, double scale) {
    FileReader file(filename);
    std::string_view line;

    std::vector<dvec3> vertices;
    imesh3 mesh;
    std::cmatch m;

    std::regex element_regex(R"(element (vertex|face) (\d+)\s*)");
    while (true) {
		line = file.readline();
        if (line == "")
            THROW(runtime_error, "bad ply file header");
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

    std::regex vertex_regex(format(R"(^(%s) (%s) (%s)(\s+|$))", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (vertices.size() < vertices.capacity()) {
		line = file.readline();
        if (line == "")
            THROW(runtime_error, "unexpected end of ply file");
        if (!search(line, vertex_regex, /*out*/m))
            THROW(runtime_error, "expected vertex in ply file: [%s]", line);
        double a = parse<double>(m[1]);
        double b = parse<double>(m[3]);
		double c = parse<double>(m[5]);
        vertices.push_back(dvec3{a, b, c} * scale);
    }

    std::regex face_regex(R"(3 (\d+) (\d+) (\d+)\s*)");
    while (mesh.size() < mesh.capacity()) {
        line = file.readline();
        if (line == "")
            THROW(runtime_error, "unexpected end of ply file");
		if (!match(line, face_regex, /*out*/m))
            THROW(runtime_error, "expected face in ply file");
        int a = parse<int>(m[1]);
        int b = parse<int>(m[3]);
		int c = parse<int>(m[5]);
		ivec3 va = vconvert(vertices[a], ivec3);
		ivec3 vb = vconvert(vertices[b], ivec3);
		ivec3 vc = vconvert(vertices[c], ivec3);
        mesh.emplace_back(va, vb, vc);
    }

    return mesh;
}
