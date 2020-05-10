#include <geom/segment.h>
#include <geom/tesselate.h>

static bool intersects(segment2 p, segment2 q) {
    double sp = signed_double_area(p.a, p.b, q.a);
    double sq = signed_double_area(p.a, p.b, q.b);
    if (sp * sq > 0) return false;
    double sa = signed_double_area(q.a, q.b, p.a);
    double sb = signed_double_area(q.a, q.b, p.b);
    return sa * sb < 0;
}

static bool intersects(const polygon2& polygon, segment2 p) {
    for (segment2 q : Edges(polygon))
        if (intersects(p, q)) return true;
    return false;
}

void tesselate(polygon2 polygon, mesh2& mesh) {
    uint b = 0;
    double a1 = signed_double_area(polygon);
    while (polygon.size() > 3) {
        b += 1;
        if (b >= polygon.size()) b = 0;

        uint a = (b > 0) ? b - 1 : (polygon.size() - 1);
        uint c = (b < polygon.size() - 1) ? b + 1 : 0;
        segment2 e(polygon[a], polygon[c]);

        double a2 = signed_double_area(e.a, polygon[b], e.b);
        if (abs(abs(a2) + abs(a1 - a2) - abs(a1)) >= 1e-9 || intersects(polygon, e)) continue;

        mesh.emplace_back(e.a, polygon[b], e.b);
        a1 -= a2;
        polygon.erase(polygon.begin() + b);
    }
    mesh.emplace_back(polygon[0], polygon[1], polygon[2]);
}
