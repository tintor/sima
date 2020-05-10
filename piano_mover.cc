#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>
#include <geom/classify.h>
#include <geom/pose.h>
#include <geom/properties.h>
#include <view/font.h>
#include <view/glm.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>
#include <view/window.h>

constexpr int Width = 1000, Height = 1000;

constexpr double LinearSpeed = 100;
constexpr double AngularSpeed = M_PI / 10;

struct Info {
    pose2 prev;
    double dist;  // direct distance (not total dist)
};

struct Model {
    std::default_random_engine rnd;

    // input
    aabb2 space;
    vector<double2> obstacles;  // points
    vector<double2> piano;      // polygon (for zero pos and zero heading)
    pose2 start;
    pose2 end;

    // output
    vector<pose2> solution;
    vector<double> solution_dist;
    size_t step;
    double step_dist = 0;
    bool paused = false;

    unordered_map<pose2, Info, hash_t<pose2>> graph;

    Model() : rnd(std::random_device()()) { init(); }

    void init() {
        space.min = {0, 0};
        space.max = {Width, Height};

        piano.resize(4);
        piano[0] = {50, 10};
        piano[1] = {-50, 10};
        piano[2] = {-50, -10};
        piano[3] = {50, -10};

        start = pose2(double2{100, 100}, 0);
        end = pose2(double2{900, 900}, M_PI / 2);

        auto piano_s = transform(piano, start);
        auto piano_e = transform(piano, end);
        obstacles.clear();
        for (size_t i = 0; i < 100; i++)
            while (true) {
                double2 p = uniform2(rnd, space);
                if (!PointInPolygon(p, piano_s) && !PointInPolygon(p, piano_e)) {
                    obstacles.push_back(p);
                    break;
                }
            }

        solution.clear();
        step = 0;

        graph.clear();
        graph[start] = {start, 0};
    }

    static vector<double2> transform(const vector<double2>& poly, pose2 pose) {
        vector<double2> result(poly.size());
        for (size_t i = 0; i < poly.size(); i++) result[i] = pose.apply(poly[i]);
        return result;
    }

    pose2 sample_valid_pose() {
        while (true) {
            pose2 p;
            p.position = uniform2(rnd, space);
            p.orientation = uniform(rnd, -M_PI, M_PI);
            if (!collision(p)) return p;
        }
    }

    bool collision(pose2 pose) {
        auto piano_w = transform(piano, pose);
        for (double2 obs : obstacles)
            if (PointInPolygon(obs, piano_w)) return true;
        return false;
    }

    // time it takes to move from src to dest
    static double distance(pose2 src, pose2 dest) {
        return max(length(src.position - dest.position) / LinearSpeed, angle(src, dest) / AngularSpeed);
    }

    bool is_valid_move(pose2 src, pose2 dest) {
        constexpr int N = 1000;
        for (int i = 0; i <= N; i++)
            if (collision(interpolate(src, dest, double(i) / N))) return false;
        return true;
    }

    pose2 closest_pose_in_graph(pose2 p) const {
        double dist_min = std::numeric_limits<double>::max();
        const pose2* pose_min;

        for (auto it = graph.begin(); it != graph.end(); it++) {
            double dist = distance(it->first, p);
            if (dist < dist_min) {
                dist_min = dist;
                pose_min = &it->first;
            }
        }
        return *pose_min;
    }

    // move from src to dest until collision
    pose2 truncate_pose(pose2 src, pose2 dest) {
        double d = distance(src, dest);
        /*if (d > 1) {
                dest = interpolate(src, dest, 1 / d);
                d = 1;
        }*/

        constexpr double dt = 0.0001;
        size_t steps = ceil(d / dt);
        pose2 valid = src;
        for (size_t i = 1; i <= steps; i++) {
            pose2 mid = interpolate(src, dest, double(i) / steps);
            if (collision(mid)) break;
            valid = mid;
        }
        return valid;
    }

    double start_dist(pose2 a) const {
        size_t z = 0;
        double dist = 0;
        while (a != start) {
            const auto& i = graph.at(a);
            dist += i.dist;
            a = i.prev;
            if (++z >= 1000000) THROW(runtime_error, "graph loop");
        }
        return dist;
    }

    void improve_graph(pose2 p) {
        constexpr double MaxDist = 1;
        // TODO only fetch poses within MaxDist of p
        double dist_s_p = start_dist(p);
        for (auto it = graph.begin(); it != graph.end(); it++) {
            pose2 m = it->first;
            auto& i = it->second;

            if (m == p || m == start) continue;
            double dist_p_m = distance(p, m);
            if (dist_p_m < MaxDist && i.dist > 0 && dist_s_p + dist_p_m < start_dist(i.prev) + i.dist &&
                is_valid_move(p, m)) {
                i.prev = p;
                i.dist = dist_p_m;
            }
        }
    }

    void solve(size_t steps) {
        for (size_t i = 0; i < steps; i++) {
            pose2 p = sample_valid_pose();
            pose2 m = closest_pose_in_graph(p);
            p = truncate_pose(m, p);
            if (p != m) {
                // reduce parent
                double dist_m_p = distance(m, p);
                /*while (m != start) {
                        pose2 e = graph[m].prev;
                        if (!is_valid_move(e, p))
                                break;
                        double dist_e_p = distance(e, p);
                        if (dist_e_p >= graph[m].dist + dist_m_p)
                                break;
                        m = e;
                        dist_m_p = dist_e_p;
                }*/

                graph[p] = {m, dist_m_p};
                improve_graph(p);

                double dist_p_e = distance(p, end);
                if (dist_p_e <= 1 && is_valid_move(p, end)) {
                    // TOOD existing solution might be better than PE
                    graph[end] = {p, dist_p_e};
                    step = 0;
                    solution.clear();
                    p = end;
                    while (p != start) {
                        solution.push_back(p);
                        p = graph[p].prev;
                    }
                    solution.push_back(start);
                    std::reverse(solution.begin(), solution.end());
                    solution_dist.clear();
                    solution_dist.push_back(0);
                    for (size_t i = 1; i < solution.size(); i++)
                        solution_dist.push_back(solution_dist.back() + graph[solution[i]].dist);
                    break;
                }
            }
        }
    }
} model;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    const char* key_name = glfwGetKeyName(key, 0);
    // print("key_callback [%s] key:%s scancode:%s action:%s mods:%s\n", key_name, key, scancode, action, mods);

    if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT)
        glfwSetWindowShouldClose(window, GL_TRUE);

    if (action == GLFW_PRESS && key == GLFW_KEY_LEFT)
        if (model.step > 0) model.step -= 1;
    if (action == GLFW_PRESS && key == GLFW_KEY_LEFT)
        if (model.step + 1 < model.solution.size()) model.step += 1;

    if (action == GLFW_PRESS && key == GLFW_KEY_ENTER) model.init();

    if (action == GLFW_PRESS && key == GLFW_KEY_SPACE) model.paused = !model.paused;

    if (action == GLFW_PRESS && key == GLFW_KEY_1) model.solve(1);
    if (action == GLFW_PRESS && key == GLFW_KEY_2) model.solve(10);
    if (action == GLFW_PRESS && key == GLFW_KEY_3) model.solve(100);
    if (action == GLFW_PRESS && key == GLFW_KEY_4) model.solve(1000);
    if (action == GLFW_PRESS && key == GLFW_KEY_5) model.solve(10000);
}

void render_piano(uint64_t color, const mat3& ortho, pose2 pose, Uniform_mat3& transformUniform,
                  VertexBuffer_vec2_rgba& buffer) {
    mat3 transform(ortho);
    transform = glm::translate(transform, vec2(pose.position.x, pose.position.y));
    transform = glm::rotate(transform, float(-pose.orientation));
    transformUniform = transform;

    for (const auto& p : model.piano) buffer.add(p, color);
    buffer.draw(GL_LINE_LOOP);
}

int main(int argc, char** argv) {
    InitSegvHandler();

    auto window = CreateWindow({.width = Width, .height = Height, .resizeable = false});
    glClearColor(0.0, 0.0, 0.0, 1.0);

    glfwSetFramebufferSizeCallback(window,
                                   [](GLFWwindow* window, int width, int height) { glViewport(0, 0, width, height); });

    glfwSetKeyCallback(window, key_callback);

    Shader shader(R"END(
		#version 330 core
		layout (location = 0) in vec2 pos;
		layout (location = 1) in vec4 color;
		out vec4 color_f;
		uniform mat3 transform;

		void main() {
		    vec3 p = transform * vec3(pos, 1.0);
		    gl_Position = vec4(p.x, p.y, 0.0, 1.0);
			color_f = color;
		}

		#version 330 core
		out vec4 color;
		in vec4 color_f;

		void main() {
		    color = color_f;
		}
	)END");
    UNIFORM(mat3, transform);

    FontRenderer renderer(800, 600);
    Font monaco("/System/Library/Fonts/Monaco.dfont", 48, &renderer);

    mat3 ortho;
    ortho[0][0] = 2.0f / Width;
    ortho[1][1] = 2.0f / Height;
    ortho[2][2] = 1.0f;
    ortho[2][0] = -1.0f;
    ortho[2][1] = -1.0f;

    VertexBuffer_vec2_rgba buffer(2000);
    glPointSize(2);

    RunEventLoop(window, [&]() {
        if (!model.paused) {
            if (model.solution.size() > 0) {
                model.step_dist += 0.025;
            } else {
                model.step = 0;
                model.step_dist = 0;
            }
        }

        glClear(GL_COLOR_BUFFER_BIT);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        monaco.render("hello!", 5, 5, 0.3, "7FE030");
        glDisable(GL_BLEND);

        glUseProgram(shader);
        buffer.bind();

        // render obstacles
        mat3 transform(ortho);
        transformUniform = transform;
        for (const auto& obs : model.obstacles) buffer.add(obs, 0xFFFFFFFF);
        buffer.draw(GL_POINTS);

        // render graph
        for (auto it = model.graph.begin(); it != model.graph.end(); it++) {
            buffer.add(it->first.position, 0xFF0000FF);
            buffer.add(it->second.prev.position, 0xFF888888);
            if (buffer.data.size() == 2000) buffer.draw(GL_LINES);
        }
        buffer.draw(GL_LINES);

        // render piano
        render_piano(0xFF0000FF, ortho, model.end, transformUniform, buffer);
        render_piano(0xFF00FF00, ortho, model.start, transformUniform, buffer);
        if (model.solution.size() >= 2) {
            auto& s = model.step;
            const auto& sol = model.solution;
            const auto& sd = model.solution_dist;

            if (model.step_dist >= sd.back()) {
                model.step_dist = 0;
                model.step = 0;
            }
            while (s + 1 < sd.size() && sd[s + 1] <= model.step_dist) s += 1;
            double k = (model.step_dist - sd[s]) / (sd[s + 1] - sd[s]);
            auto pose = interpolate(sol[s], sol[s + 1], k);
            render_piano(0xFFFFFFFF, ortho, pose, transformUniform, buffer);
        }
    });
    return 0;
}
