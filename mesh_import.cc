#include "mesh_import.h"
#include "util.h"
#include "file.h"
#include <regex>

constexpr const char* FLOAT_REGEX = R"(-?\d+\.\d+(e[+-]\d+)?)";

Mesh3d load_stl(const std::string& filename) {
    FileReader file(filename.c_str());
    std::string_view line;

    Mesh3d mesh;
    std::cmatch match;

    if (!file.getline(line) || line != "solid Default")
        throw std::runtime_error("expected [solid Default]");

    std::regex facet_regex(R"(^\s*facet\s+)");
    std::regex vertex_regex(tfm::format(R"(^\s*vertex (%s) (%s) (%s)\s*$)", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (true) {
        if (!file.getline(line))
            throw std::runtime_error("expected [facet ...] or [endsolid Default] insted of end");
        if (line == "endsolid Default")
            break;
        if (!std::regex_search(line.begin(), line.end(), facet_regex))
            throw std::runtime_error(tfm::format("expected [facet ...] or [endsolid Default] instead of [%s]", line));
        if (!file.getline(line) || line != "    outer loop")
            throw std::runtime_error("expected [outer loop]");

        dvec3 v[3];
        FOR(i, 3) {
            if (!file.getline(line) || !std::regex_match(line.begin(), line.end(), match, vertex_regex))
                throw std::runtime_error(tfm::format("expected [vertex ...] instead of [%s]", line));

            from_chars(match[1].first, match[1].second, v[i].x);
            from_chars(match[3].first, match[3].second, v[i].y);
            from_chars(match[5].first, match[5].second, v[i].z);
        }
        mesh.emplace_back(v[0], v[1], v[2]);

        if (!file.getline(line) || line != "    endloop")
            throw std::runtime_error("expected [endloop]");
        if (!file.getline(line) || line != "  endfacet")
            throw std::runtime_error("expected [endfacet]");
    }

    return mesh;
}

Mesh3d load_ply(const std::string& filename) {
    FileReader file(filename.c_str());
    std::string_view line;

    std::vector<dvec3> vertices;
    Mesh3d mesh;
    std::cmatch match;

    std::regex element_regex(R"(element (vertex|face) (\d+)\s*)");
    while (true) {
        if (!file.getline(line))
            throw std::runtime_error("bad ply file header");
        if (line == "end_header")
            break;
        if (std::regex_match(line.begin(), line.end(), /*out*/match, element_regex)) {
            int s;
            from_chars(match[2].first, match[2].second, s);
            if (match[1].str() == "vertex")
                vertices.reserve(s);
            else
                mesh.reserve(s);
        }
    }

    std::regex vertex_regex(tfm::format(R"(^(%s) (%s) (%s)(\s+|$))", FLOAT_REGEX, FLOAT_REGEX, FLOAT_REGEX));
    while (vertices.size() < vertices.capacity()) {
        if (!file.getline(line))
            throw std::runtime_error("unexpected end of ply file");
        if (!std::regex_search(line.begin(), line.end(), /*out*/match, vertex_regex))
            throw std::runtime_error(tfm::format("expected vertex in ply file: [%s]", line));
        double a, b, c;
        from_chars(match[1].first, match[1].second, a);
        from_chars(match[3].first, match[3].second, b);
        from_chars(match[5].first, match[5].second, c);
        vertices.emplace_back(a, b, c);
    }

    std::regex face_regex(R"(3 (\d+) (\d+) (\d+)\s*)");
    while (mesh.size() < mesh.capacity()) {
        if (!file.getline(line) || !std::regex_match(line.begin(), line.end(), /*out*/match, face_regex))
            throw std::runtime_error("expected face in ply file");
        int a, b, c;
        from_chars(match[1].first, match[1].second, a);
        from_chars(match[3].first, match[3].second, b);
        from_chars(match[5].first, match[5].second, c);
        mesh.emplace_back(vertices[a], vertices[b], vertices[c]);
    }

    return mesh;
}
