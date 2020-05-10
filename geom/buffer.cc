#include <core/exception.h>
#include <geom/buffer.h>

static void GenerateCircularArc(vector<double2>& arc, double2 center, double start_angle, double end_angle,
                                double radius, int count, bool end) {
    ASSERT_ALWAYS(start_angle < end_angle, "%s < %s", start_angle, end_angle);
    ASSERT_ALWAYS(radius >= 0);
    ASSERT_ALWAYS(count >= 0);
    double delta = (end_angle - start_angle) / count;
    arc.push_back(center + double2{cos(start_angle), sin(start_angle)} * radius);
    for (int i = 1; i < count; i++) {
        double angle = start_angle + i * delta;
        arc.push_back(center + double2{cos(angle), sin(angle)} * radius);
    }
    if (end) arc.push_back(center + double2{cos(end_angle), sin(end_angle)} * radius);
}

// returns angle in range [-PI, PI]
static double angle(double2 v) { return std::atan2(v.y, v.x); }

vector<double2> ComputeBuffer(cspan<double2> convex, double radius, int circle_vertices) {
    ASSERT_ALWAYS(radius > 0);
    ASSERT_ALWAYS(circle_vertices >= 4);
    ASSERT_ALWAYS(convex.size() > 0);

    vector<double2> result;
    if (convex.size() == 1) {
        GenerateCircularArc(result, convex[0], 0, 2 * PI, radius, circle_vertices, false);
        return result;
    }

    // Assume the polygon is convex and counter clockwise
    const double k = circle_vertices / (2 * PI);
    const int n = convex.size();
    for (int i = 0; i < n; i++) {
        double2 prev = convex[(i == 0) ? n - 1 : (i - 1)];
        double2 curr = convex[i];
        double2 next = convex[(i == n - 1) ? 0 : (i + 1)];

        double start_angle = angle(prev - curr);
        double end_angle = angle(next - curr);
        end_angle -= PI / 2;
        start_angle += PI / 2;
        if (end_angle < start_angle) end_angle += 2 * PI;

        int points = std::round((end_angle - start_angle) * k);
        GenerateCircularArc(result, curr, start_angle, end_angle, radius, points, true);
    }
    return result;
}
