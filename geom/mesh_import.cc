#include <geom/mesh_import.h>
#include <core/util.h>
#include <core/file.h>
#include <core/range.h>
#include <core/exception.h>
#include <core/string_util.h>

constexpr const char* FLOAT_REGEX = R"(-?\d+\.\d+(e[+-]\d+)?)";

mesh3 load_stl(string_view filename) {
    FileReader file(filename);
    string_view line;

    mesh3 mesh;
    std::cmatch m;

    if (file.readline() != "solid Default\n")
        THROW(runtime_error, "expected [solid Default]");

    regex facet_regex(R"(^\s*facet\s+)");
    regex vertex_regex(format(R"(^\s*vertex (%s) (%s) (%s)\s*$)", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (true) {
		line = file.readline();
        if (line == "endsolid Default\n")
            break;
        if (!search(line, facet_regex))
            THROW(runtime_error, "expected [facet ...] or [endsolid Default] instead of [%s]", line);
        if (file.readline() != "    outer loop\n")
            THROW(runtime_error, "expected [outer loop]");

        double3 v[3];
        for (auto i : range(3)) {
			line = file.readline();
            if (!match(line, vertex_regex, /*out*/m))
                THROW(runtime_error, "expected [vertex ...] instead of [%s]", line);
            v[i].x = parse<double>(m[1]);
            v[i].y = parse<double>(m[3]);
            v[i].z = parse<double>(m[5]);
        }
        mesh.emplace_back(v[0], v[1], v[2]);

        if (file.readline() != "    endloop\n")
            THROW(runtime_error, "expected [endloop]");
        if (file.readline() != "  endfacet\n")
            THROW(runtime_error, "expected [endfacet]");
    }

    return mesh;
}

mesh3 load_ply(string_view filename) {
    FileReader file(filename);
    string_view line;

    vector<double3> vertices;
    mesh3 mesh;
    std::cmatch m;

    regex element_regex(R"(element (vertex|face) (\d+)\s*)");
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

    regex vertex_regex(format(R"(^(%s) (%s) (%s)(\s+|$))", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (vertices.size() < vertices.capacity()) {
		line = file.readline();
        if (line == "")
            THROW(runtime_error, "unexpected end of ply file");
        if (!search(line, vertex_regex, /*out*/m))
            THROW(runtime_error, "expected vertex in ply file: [%s]", line);
        double a = parse<double>(m[1]);
        double b = parse<double>(m[3]);
		double c = parse<double>(m[5]);
        vertices.push_back(double3{a, b, c});
    }

    regex face_regex(R"(3 (\d+) (\d+) (\d+)\s*)");
    while (mesh.size() < mesh.capacity()) {
        line = file.readline();
        if (line == "")
            THROW(runtime_error, "unexpected end of ply file");
		if (!match(line, face_regex, /*out*/m))
            THROW(runtime_error, "expected face in ply file");
        int a = parse<int>(m[1]);
        int b = parse<int>(m[3]);
		int c = parse<int>(m[5]);
        mesh.emplace_back(vertices[a], vertices[b], vertices[c]);
    }

    return mesh;
}
