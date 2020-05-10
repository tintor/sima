#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>
#include <view/font.h>
#include <view/glm.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>
#include <view/window.h>

#include <random>
#include <variant>

// Board and judge
// ===============

enum class Figure : char { None, Dome, Player1, Player2 };

struct Cell {
    char level = 0;
    Figure figure = Figure::None;
};

struct Coord {
    int x, y;
};

bool operator==(Coord a, Coord b) { return a.x == b.x && a.y == b.y; }
bool operator!=(Coord a, Coord b) { return !(a == b); }
bool Nearby(Coord src, Coord dest) { return abs(src.x - dest.x) <= 1 && abs(src.y - dest.y) <= 1; }
bool IsValid(Coord a) { return a.x >= 0 && a.x < 5 && a.y >= 0 && a.y < 5; }

struct Board {
    bool setup = true;
    Figure player = Figure::Player1;
    optional<Coord> moved;
    bool built = false;

    std::array<std::array<Cell, 5>, 5> cell;

    const Cell& operator()(Coord c) const { return cell[c.x][c.y]; }
    Cell& operator()(Coord c) { return cell[c.x][c.y]; }
};

Board g_board;

Figure Other(Figure player) { return (player == Figure::Player1) ? Figure::Player2 : Figure::Player1; }

vector<Coord> All() {
    vector<Coord> out;
    out.reserve(25);
    for (int x = 0; x < 5; x++)
        for (int y = 0; y < 5; y++) out.push_back({x, y});
    return out;
}

const vector<Coord> kAll = All();

template <typename Fn>
int Count(const Board& board, const Fn& fn) {
    int c = 0;
    for (Coord e : kAll)
        if (fn(board(e))) c += 1;
    return c;
}

#define L(A) [&](const auto& e) { return A; }

struct NextAction {};
struct PlaceAction {
    Coord dest;
};
struct MoveAction {
    Coord src, dest;
};
struct BuildAction {
    Coord dest;
    bool dome;
};

using Action = std::variant<NextAction, PlaceAction, MoveAction, BuildAction>;

optional<string_view> Next(Board& board) {
    if (board.setup) {
        if (Count(board, L(e.figure == board.player)) != 2) return "need to place worker";
    } else {
        if (!board.moved) return "need to move";
        if (!board.built) return "need to build";
    }

    board.player = Other(board.player);
    board.moved = std::nullopt;
    board.built = false;
    if (board.setup && Count(board, L(e.figure != Figure::None)) == 4) board.setup = false;
    return nullopt;
}

optional<string_view> Place(Board& board, Coord dest) {
    if (!IsValid(dest)) return "invalid coord";
    if (!board.setup) return "can't place after setup is complete";
    if (board(dest).figure != Figure::None) return "occupied";
    if (Count(board, L(e.figure == board.player)) == 2) return "can't place anymore";

    board(dest).figure = board.player;
    return nullopt;
}

optional<string_view> Move(Board& board, Coord src, Coord dest) {
    if (!IsValid(src) || !IsValid(dest)) return "invalid coord";
    if (board.setup) return "can't move during setup";
    if (board.moved) return "moved already";
    if (board(src).figure != board.player) return "player doesn't have figure at src";
    if (board(dest).figure != Figure::None) return "dest isn't empty";
    if (!Nearby(src, dest)) return "src and dest aren't nearby";
    if (board(dest).level - board(src).level > 1) return "dest is too high";

    board(dest).figure = board.player;
    board(src).figure = Figure::None;
    board.moved = dest;
    return nullopt;
}

optional<string_view> Build(Board& board, Coord dest, bool dome) {
    if (!IsValid(dest)) return "invalid coord";
    if (board.setup) return "can't build during setup";
    if (!board.moved) return "need to move";
    if (board.built) return "already built";
    if (board(dest).figure != Figure::None) return "can only build on empty space";
    if (dome && board(dest).level != 3) return "dome can only be built on level 3";
    if (!dome && board(dest).level == 3) return "floor can only be built on levels 0, 1 and 2";
    if (!Nearby(*board.moved, dest)) return "can only build near moved figure";

    if (dome) {
        board(dest).figure = Figure::Dome;
    } else {
        board(dest).level += 1;
    }
    board.built = true;
    return nullopt;
}

template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...)->overloaded<Ts...>;  // not needed as of C++20

optional<string_view> Execute(Board& board, const Action& action) {
    return std::visit(
        overloaded{[&](NextAction a) { return Next(board); }, [&](PlaceAction a) { return Place(board, a.dest); },
                   [&](MoveAction a) { return Move(board, a.src, a.dest); },
                   [&](BuildAction a) { return Build(board, a.dest, a.dome); }},
        action);
}

void Print(const Action& action) {
    std::visit(overloaded{[&](NextAction a) { print("next"); },
                          [&](PlaceAction a) { print("place:%s%s", a.dest.x, a.dest.y); },
                          [&](MoveAction a) { print("move:%s%s:%s%s", a.src.x, a.src.y, a.dest.x, a.dest.y); },
                          [&](BuildAction a) { print("build:%s%s:%s", a.dest.x, a.dest.y, a.dome ? 'D' : 'T'); }},
               action);
}

bool IsMoveBlocked(const Board& board) {
    for (Coord a : kAll)
        if (board(a).figure == board.player) {
            Board board2 = board;
            for (Coord b : kAll) {
                if (Move(board2, a, b) == nullopt) return false;
            }
        }
    return true;
}

bool IsBuildBlocked(const Board& board) {
    if (!board.moved) return false;
    Board board2 = board;
    for (Coord b : kAll)
        if (Build(board2, b, false) == nullopt || Build(board2, b, true) == nullopt) return false;
    return true;
}

bool OnThirdLevel(const Board& board) { return Count(board, L(e.figure == board.player && e.level == 3)) > 0; }

Figure Winner(const Board& board) {
    if (board.setup) return Figure::None;

    if (!board.moved && IsMoveBlocked(board)) return Other(board.player);
    if (!board.built && IsBuildBlocked(board)) return Other(board.player);
    if (OnThirdLevel(board)) return board.player;
    return Figure::None;
}

// Computer interface
// ==================

std::random_device rd;
std::mt19937 g_random(rd());

int RandomInt(int count) { return std::uniform_int_distribution<int>(0, count - 1)(g_random); }

Coord RandomCoord() {
    std::uniform_int_distribution<int> dis(0, 4);
    return Coord{dis(g_random), dis(g_random)};
}

vector<Coord> MyFigures(const Board& board) {
    vector<Coord> out;
    for (Coord e : kAll)
        if (board(e).figure == board.player) out.push_back(e);
    return out;
}

Coord MyRandomFigure(const Board& board) {
    Coord out;
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    int count = 0;
    for (Coord e : kAll)
        if (board(e).figure == board.player && dis(g_random) <= 1.0 / ++count) out = e;
    return out;
}

Action RandomAction(const Board& board) {
    if (board.setup) {
        switch (RandomInt(2)) {
            case 0:
                return NextAction{};
            case 1:
                return PlaceAction{.dest = RandomCoord()};
        }
    }

    switch (RandomInt(3)) {
        case 0:
            return NextAction{};
        case 1:
            return MoveAction{.src = MyRandomFigure(board), .dest = RandomCoord()};
        case 2:
            return BuildAction{.dest = RandomCoord(), .dome = bool(RandomInt(2))};
    }
    throw std::runtime_error("unreachable");
}

bool IsValid(const Board& board, const Action& action) {
    Board my_board = board;
    return Execute(my_board, action) == nullopt;
}

template <typename Visitor>
bool Visit(const Board& board, const Action& action, const Visitor& visit) {
    Board my_board = board;
    if (Execute(my_board, action) != nullopt) return true;
    return visit(my_board, action);
}

#define VISIT(A) \
    if (!Visit(board, A, visit)) return false;

template <typename Visitor>
bool AllValidActions(const Board& board, const Visitor& visit) {
    VISIT(NextAction{});
    if (board.setup) {
        for (Coord e : kAll) VISIT(PlaceAction{e});
        return true;
    }
    for (Coord e : kAll) {
        if (board(e).figure == board.player) {
            for (Coord d : kAll)
                if (d != e) VISIT((MoveAction{e, d}));
        }
    }
    for (bool dome : {false, true})
        for (Coord e : kAll) VISIT((BuildAction{e, dome}));
    return true;
}

template <typename Visitor>
bool AllValidActionSequences(const Board& board, vector<Action>& temp, const Visitor& visit) {
    return AllValidActions(board, [&](const Board& new_board, const Action& action) {
        temp.push_back(action);
        auto winner = Winner(new_board);
        if (std::holds_alternative<NextAction>(action) || winner != Figure::None) {
            if (!visit(temp, winner)) return false;
        } else {
            if (!AllValidActionSequences(new_board, temp, visit)) return false;
        }
        temp.pop_back();
        return true;
    });
}

void Print(const vector<Action>& actions) {
    for (const Action& action : actions) {
        Print(action);
        print(" ");
    }
    print("\n");
}

Action AutoGreedy(const Board& board) {
    vector<Action> temp;
    Action choice;
    size_t count = 0;
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    AllValidActionSequences(board, temp, [&](const vector<Action>& actions, Figure winner) {
        // Print(actions);
        if (winner == board.player) {
            choice = actions[0];
            return false;
        }
        if (winner == Figure::None && dis(g_random) <= 1.0 / ++count) choice = actions[0];
        return true;
    });
    return choice;
}

Action AutoRandom(const Board& board) {
    while (true) {
        Action action = RandomAction(board);
        if (IsValid(board, action)) return action;
    }
}

// Human interface
// ===============

constexpr int Width = 1000, Height = 1000;

vector<Board> g_history;
Board g_board_copy;
optional<Coord> g_selected;

string_view PlayerName(Figure player) { return (player == Figure::Player1) ? "yellow"sv : "red"sv; }

void Play(optional<string_view> status) {
    if (status) {
        print("%s\n", *status);
    } else {
        g_history.push_back(g_board_copy);
        g_selected = nullopt;
        auto w = Winner(g_board);
        if (w != Figure::None) {
            print("Player %s wins!\n", PlayerName(w));
        }
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
        glfwSetWindowShouldClose(window, GL_TRUE);
        return;
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_ENTER) {
        g_board_copy = g_board;
        Play(Next(g_board));
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_U) {
        if (g_history.size() > 0) {
            g_board = g_history.back();
            g_history.pop_back();
        }
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_1) {
        g_board_copy = g_board;
        Action action = AutoRandom(g_board);
        Play(Execute(g_board, action));
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_2) {
        g_board_copy = g_board;
        print("auto greedy begin\n");
        Action action = AutoGreedy(g_board);
        print("auto greedy end\n");
        Play(Execute(g_board, action));
    }
}

void OnClick(Coord dest, bool left, bool shift) {
    g_board_copy = g_board;
    if (!left) {
        Play(Build(g_board, dest, shift));
    } else if (g_board.setup) {
        Play(Place(g_board, dest));
    } else if (g_board(dest).figure == g_board.player) {
        g_selected = dest;
    } else if (g_selected) {
        Play(Move(g_board, *g_selected, dest));
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS) {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        int mouse_x = floor(xpos - 100) / 160;
        int mouse_y = floor(Height - ypos - 100) / 160;
        if (mouse_x >= 0 && mouse_x < 5 && mouse_y >= 0 && mouse_y < 5)
            OnClick({mouse_x, mouse_y}, button == GLFW_MOUSE_BUTTON_LEFT, mods == GLFW_MOD_SHIFT);
    }
}

void scroll_callback(GLFWwindow* window, double x, double y) {}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) { glViewport(0, 0, width, height); }

struct View : public VertexBuffer_vec2_rgba {
    Shader shader;
    Uniform_mat4 transform;

    FontRenderer font_renderer;
    Font mono;
    glm::mat4 ortho;

    View()
        : VertexBuffer_vec2_rgba(25),
          shader(R"END(
    		#version 330 core
    		layout (location = 0) in vec2 pos;
    		layout (location = 1) in vec4 rgba;
    		out vec4 pixel_rgba;
    		uniform mat4 transform;

            void main() {
    		    vec4 p = transform * vec4(pos, 0.0, 1.0);
		        gl_Position = p;
                pixel_rgba = rgba;
	    	}

		    #version 330 core
    		in vec4 pixel_rgba;
    		out vec4 color;

    		void main() {
    		    color = pixel_rgba;
    		}
        	)END"),
          font_renderer(1000, 1000),
          mono("JetBrainsMono-Medium.ttf", 48, &font_renderer),
          transform("transform") {}
};

void Render(const Board& board, View& view) {
    glUseProgram(view.shader);
    view.bind();
    view.transform = view.ortho;

    // Frame
    for (int i = 0; i <= 5; i++) {
        const uint64_t color = 0xFF808080;
        view.add({100, 100 + i * 160}, color);
        view.add({100 + 800, 100 + i * 160}, color);

        view.add({100 + i * 160, 100}, color);
        view.add({100 + i * 160, 100 + 800}, color);
    }
    view.draw(GL_LINES);

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            const double s = 160;
            const double px = 100 + x * s + s / 2;
            const double py = 100 + y * s + s / 2;
            const auto cell = board.cell[x][y];

            if (g_selected && x == g_selected->x && y == g_selected->y) {
                double a = 80;
                const uint64_t color = 0xFF00FFFF;
                view.add({px - a, py - a}, color);
                view.add({px + a, py - a}, color);
                view.add({px + a, py + a}, color);
                view.add({px - a, py + a}, color);
                view.draw(GL_LINE_LOOP);
            }

            // Towers
            if (cell.level >= 1) {
                double a = 70;
                const uint64_t color = 0xFFFFFFFF;
                view.add({px - a, py - a}, color);
                view.add({px + a, py - a}, color);
                view.add({px + a, py + a}, color);
                view.add({px - a, py + a}, color);
                view.draw(GL_LINE_LOOP);
            }
            if (cell.level >= 2) {
                const uint64_t color = 0xFFFFFFFF;
                double a = 60;
                view.add({px - a, py - a}, color);
                view.add({px + a, py - a}, color);
                view.add({px + a, py + a}, color);
                view.add({px - a, py + a}, color);
                view.draw(GL_LINE_LOOP);
            }
            if (cell.level >= 3) {
                const uint64_t color = 0xFFFFFFFF;
                double a = 50;
                view.add({px - a, py - a}, color);
                view.add({px + a, py - a}, color);
                view.add({px + a, py + a}, color);
                view.add({px - a, py + a}, color);
                view.draw(GL_LINE_LOOP);
            }

            // Domes
            if (cell.figure == Figure::Dome) {
                const uint64_t color = 0xFFFF0000;
                double a = 40;
                for (int i = 0; i < 16; i++) {
                    double k = i * 2 * M_PI / 16;
                    view.add({px + cos(k) * a, py + sin(k) * a}, color);
                }
                view.draw(GL_TRIANGLE_FAN);
            }

            // Player pieces
            if (cell.figure == Figure::Player1 || cell.figure == Figure::Player2) {
                const uint64_t color = (cell.figure == Figure::Player1) ? 0xFF00FFFF : 0xFF0000FF;
                double a = 40;
                view.add({px - a, py - a}, color);
                view.add({px + a, py - a}, color);
                view.add({px + a, py + a}, color);
                view.add({px - a, py + a}, color);
                view.draw(GL_TRIANGLE_FAN);
            }
        }
    }

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    view.mono.moveTo(10, 980);
    view.mono.m_scale = 13.0 / 48;
    view.mono.m_color = Color("FFFFFF");
    view.mono.render(format("Player %s", PlayerName(board.player)));
    if (board.moved) view.mono.render(" Moved");
    if (board.built) view.mono.render(" Built");
    glDisable(GL_BLEND);
}

// u - undo
// place worker - left click
// select worker - left click
// move worker - left click
// build tower - right click
// build dome - shift + right click
// Enter - done

using Policy = std::function<Action(const Board&)>;

Figure Battle(const Policy& policy_a, const Policy& policy_b) {
    Board board;
    while (true) {
        auto w = Winner(board);
        if (w != Figure::None) return w;
        const Policy& policy = (board.player == Figure::Player1) ? policy_a : policy_b;
        auto s = Execute(board, policy(board));
        if (s != nullopt) {
            print("faul\n");
            return Other(board.player);
        }
    }
}

void AutoBattle(int count, string_view name_a, const Policy& policy_a, string_view name_b, const Policy& policy_b) {
    int wins_a = 0;
    for (int i = 0; i < count; i++)
        if (Battle(policy_a, policy_b) == Figure::Player1) wins_a += 1;
    print("%s %s : %s %s\n", name_a, wins_a, count - wins_a, name_b);
}

int main(int argc, char** argv) {
    InitSegvHandler();
    AutoBattle(100, "random", AutoRandom, "greedy", AutoGreedy);
    AutoBattle(100, "greedy", AutoGreedy, "random", AutoRandom);
    return 0;

    auto window = CreateWindow({.width = Width, .height = Height, .resizeable = false});
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glViewport(0, 0, Width, Height);

    glClearColor(0.0, 0.0, 0.0, 1.0);

    View view;
    view.ortho = glm::ortho(0.0, double(Width), 0.0, double(Height));

    RunEventLoop(window, [&]() {
        glClear(GL_COLOR_BUFFER_BIT);
        Render(g_board, view);
    });
    return 0;
}
