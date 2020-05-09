// load human 3d model
// render human 3d model / with phong shading
// articulate human 3d model

// render contacts
// collision between concave shapes
// rotation!
// hinges and 3 dof arm (turn off collision detection between objects that have hinge)
// draw polygons with mouse
// render velocity and ang_velocity vectors
// move and rotate objects with mouse
// tangental friction
// notification messages that disappear
// joint friction
// PID controller for the arm
// click on stop to move tip of end-effector there
// WASD to move end-effector
// control for each joint individually
// open/close end-effector

// 3D:
// balls
// convex shapes (no rotation)
// rotation
// hinges
// concave shapes

#include <core/callstack.h>
#include <core/format.h>
#include <core/util.h>
#include <geom/classify.h>
#include <geom/pose.h>
#include <geom/properties.h>
#include <geom/triangle.h>
#include <sim/integration.h>
#include <view/font.h>
#include <view/glm.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>
#include <view/window.h>

bool gSimulate = false;
bool gSimulateTick = false;
bool gGravity = true;
bool gAirDrag = false;
bool gFriction = false;

bool gSpaceKey = false;
bool gShiftKey = false;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    const char* key_name = glfwGetKeyName(key, 0);
    print("key_callback [%s] key:%s scancode:%s action:%s mods:%s\n", key_name, key, scancode, action, mods);

    if (action == GLFW_PRESS && key == GLFW_KEY_SPACE) {
	gSpaceKey = true;
    }
    if (action == GLFW_RELEASE && key == GLFW_KEY_SPACE) {
	gSpaceKey = false;
    }

    if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT) {
	glfwSetWindowShouldClose(window, GL_TRUE);
	return;
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_SPACE && mods == 0) {
	gSimulate ^= 1;
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_SPACE && mods == GLFW_MOD_SHIFT) {
	gSimulateTick ^= 1;
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_G && mods == 0) {
	gGravity ^= 1;
    }
    if (action == GLFW_PRESS && key == GLFW_KEY_A && mods == 0) {
	gAirDrag ^= 1;
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT) {
    }
}

void scroll_callback(GLFWwindow* window, double x, double y) {}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) { glViewport(0, 0, width, height); }

struct Body {
    polygon2 shape;
    vector<vec2> shapef;

    aabb2 box;
    double radius;
    double mass;

    dvec2 pos;
    dvec2 vel = vec2(0, 0);
    double ang_pos = 0;
    double ang_vel = 0;

    dvec2 saved_pos;
    dvec2 saved_vel;
    double saved_ang_pos;
    double saved_ang_vel;
};

void SaveStates(vector<Body>& bodies) {
    for (Body& body : bodies) {
	body.saved_pos = body.pos;
	body.saved_vel = body.vel;
	body.saved_ang_pos = body.ang_pos;
	body.saved_ang_vel = body.ang_vel;
    }
}

void RestoreStates(vector<Body>& bodies) {
    for (Body& body : bodies) {
	body.pos = body.saved_pos;
	body.vel = body.saved_vel;
	body.ang_pos = body.saved_ang_pos;
	body.ang_vel = body.saved_ang_vel;
    }
}

void Interact(Body& a, Body& b, polygon2& temp, vector<IContact2>& contacts) {
    dvec2 d = a.pos - b.pos;
    float r = a.radius + b.radius;
    if (glm::dot(d, d) > r * r) return;
    // translate B to A's frame
    temp.resize(b.shape.size());
    for (uint i = 0; i < b.shape.size(); i++) temp[i] = b.shape[i] + double2{d.x, d.y};

    contacts.clear();
    int c = Classify(a.shape, temp, &contacts);
    if (c > 0) return;
    if (c < 0) {
	print("penetrating objects in Interact()\n");
	exit(1);
	// TODO classify will not return contacts when there is penetration! :(
	// TODO there could be multiple contacts with different levels (or none) of penetration
	// TODO move objects slightly away to prevent penetration (while keeping total kinetic + potential energy
	// constant) (just moving objects will not work when rotations are added)
    }

    // Translate contacts from A frame to world frame
    for (IContact2& c : contacts) {
	c.sa -= double2{d.x, d.y};
	c.sb -= double2{d.x, d.y};
    }

    // Resolve all contacts with dv at once!
    double w = 0;
    for (const IContact2& contact : contacts) {
	d.x = -contact.normal.x;  // TODO check direction
	d.y = -contact.normal.y;  // TODO check direction
	double u = glm::dot(d, a.vel - b.vel);
	if (u < 0) w += -u;
    }
    if (w == 0) return;

    dvec2 wa(0, 0);
    dvec2 wb(0, 0);
    for (const IContact2& contact : contacts) {
	d.x = -contact.normal.x;  // TODO check direction
	d.y = -contact.normal.y;  // TODO check direction
	double ua = glm::dot(d, a.vel);
	double ub = glm::dot(d, b.vel);
	double u = ua - ub;
	if (u < 0) {
	    auto ma = a.mass;
	    auto mb = b.mass;
	    double va = (ua * (ma - mb) + 2 * mb * ub) / (ma + mb);
	    double vb = (ub * (mb - ma) + 2 * ma * ua) / (ma + mb);
	    wa += d * (va - ua) * -u;
	    wb += d * (vb - ub) * -u;
	}
    }

    a.vel += wa / w;
    b.vel += wb / w;
}

void SetShape(Body& body, const polygon2& poly) {
    double2 com = CenterOfMass(poly);
    body.shape.clear();
    body.shape << poly;
    for (double2& v : body.shape) v -= com;

    body.shapef.resize(body.shape.size());
    for (uint i = 0; i < body.shape.size(); i++) body.shapef[i] = vec2(body.shape[i].x, body.shape[i].y);

    double r2 = 0;
    for (double2 v : body.shape) maximize(r2, dot(v, v));
    body.radius = sqrt(r2);

    body.box = aabb2(body.shape);
}

int Classify(const Body& a, const Body& b, polygon2& temp) {
    dvec2 d = a.pos - b.pos;
    float r = a.radius + b.radius;
    if (glm::dot(d, d) > r * r) return 1;

    // translate B to A's frame
    temp.resize(b.shape.size());
    for (uint i = 0; i < b.shape.size(); i++) temp[i] = b.shape[i] + double2{d.x, d.y};
    return Classify(a.shape, temp);
}

constexpr int Width = 1200, Height = 900;

// +1 - all bodies are separate
//  0 - at least two bodies in contact (and no penetration)
// -1 - at least two bodies in penetration
int Classify(const vector<Body>& bodies, polygon2& temp) {
    int result = 1;
    for (auto i : range(bodies.size())) {
	const Body& b = bodies[i];

	// Check against walls
	int c = Sign(b.pos.y + b.box.min.y, Tolerance);
	if (c < 0) return -1;
	if (c == 0) result = 0;

	c = Sign(b.pos.x + b.box.min.x, Tolerance);
	if (c < 0) return -1;
	if (c == 0) result = 0;

	c = Sign(Width - (b.pos.x + b.box.max.x), Tolerance);
	if (c < 0) return -1;
	if (c == 0) result = 0;

	c = Sign(Height - (b.pos.y + b.box.max.y), Tolerance);
	if (c < 0) return -1;
	if (c == 0) result = 0;

	// Check against other bodies
	for (auto j : range(i + 1, bodies.size())) {
	    c = Classify(b, bodies[j], temp);
	    if (c < 0) return -1;
	    if (c == 0) result = 0;
	}
    }
    return result;
}

const dvec2 gravity(0, -100);  // mVm/s^2

void InteractWithWalls(Body& body) {
    constexpr double elasticity = 0.5;
    if (body.pos.y + body.box.min.y < 0) {
	double dip = -(body.pos.y + body.box.min.y);
	// remove penetration, but preserve total energy (if possible when d >= 0)
	double d = gravity.y * 2 * dip + body.vel.y * body.vel.y;
	body.vel.y = (d > 0) ? elasticity * sqrt(d) : 0;
	body.pos.y = -body.box.min.y;
    }
    if (body.pos.x + body.box.min.x <= 0 && body.vel.x < 0) {
	body.vel.x = -body.vel.x * elasticity;
    }
    if (body.pos.x + body.box.max.x >= Width && body.vel.x > 0) {
	body.vel.x = -body.vel.x * elasticity;
    }
    if (body.pos.y + body.box.max.y >= Height && body.vel.y > 0) {
	body.vel.y = -body.vel.y * elasticity;
    }
}

// TODO if objects end up penetrating each other, what then?
// - ignore penetration (classify needs to be able to compute contacts from penetration)
// - simulate up to contact time (resolve contacts), continue simulating
// 	 - remaining_time = dt
// 	 - start: resolve collisions
// 	 - save states of all objects
// 	 - advance all objects for remaining_time
// 	 - if any collision detected (during remaining_time):
// 	   - find collision_time TODO how?
// 	   - restore states
// 	   - advance all objects to collision_time
// 	   - remaining_time -= collision_time
// 	   - goto start
// - alternate:
// 	 - remaining_time = dt
// 	 - start: resolve collisions
// 	 - save states of all objects
// 	 - advance all objects for remaining_time
// 	 - if pedetration detected (at end time):
// 	   - find first collision_time using binary search
// 	   - restore states
// 	   - advance all objects to collision_time
// 	   - remaining_time -= collision_time
// 	   - goto start
// - alternate:
//   - simulate full dt (optionally advance until collistion time for high dv collisions)
//   - iteratively fix all penetrations (both shape v shape, and shape v wall)
//     - fix penetrations one-by-one in random order (note that fix-up may cause more penetrations)
//     - classify() will have to compute contacts for penetrations too. (complicated!)

void Advance(vector<Body>& bodies, double dt) {
    for (Body& body : bodies) {
	dvec2 acc = gGravity ? gravity : dvec2(0, 0);
	dvec4 s0 = dvec4(body.pos.x, body.pos.y, body.vel.x, body.vel.y);
	auto s1 = RungeKutta4<dvec4, double>(s0, 0, dt, [&](dvec4 s, double t) { return dvec4(dvec2(s.z, s.w), acc); });
	body.pos = dvec2(s1.x, s1.y);
	body.vel = dvec2(s1.z, s1.w);

	InteractWithWalls(body);

	if (gAirDrag) {
	    body.vel *= 0.999;
	}
    }
}

double CollisionTime(const polygon2& a, const polygon2& b, vec2 a0, vec2 a1, vec2 b0, vec2 b1) {
    // TODO
    return 0;
}

int main(int argc, char** argv) {
    InitSegvHandler();

    auto window = CreateWindow({.width = Width, .height = Height, .resizeable = false});
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);

    glClearColor(0.0, 0.5, 0.0, 1.0);

    Shader shader(R"END(
		#version 330 core
		layout (location = 0) in vec2 pos;
		uniform mat3 transform;

		void main() {
		    vec3 p = transform * vec3(pos, 1.0);
		    gl_Position = vec4(p.x, p.y, 0.0, 1.0);
		}

		#version 330 core
		out vec4 color;

		void main() {
		    color = vec4(1.0, 0.5, 0.2, 1.0);
		}
	)END");
    UNIFORM(mat3, transform);

    VertexBuffer_vec2 buffer(25);

    Font timesNewRoman("Times New Roman");
    Font arial("Arial");
    Font monaco("/System/Library/Fonts/Monaco.dfont");

    std::array<vec2, 25> v;
    for (int i = 0; i < v.size() - 1; i++) {
	auto a = 2 * PI / (v.size() - 1) * i;
	v[i].x = cos(a);
	v[i].y = sin(a);
    }
    v[24] = vec2(0, 0);
    std::array<vec2, 25> w;
    buffer.write(v);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

    vector<Body> bodies;
    bodies.push_back({.mass = 130 * 130, .radius = 130, .pos = vec2(200, 200)});
    bodies.push_back({.mass = 90 * 90, .radius = 90, .pos = vec2(400, 400)});
    bodies.push_back({.mass = 50 * 50, .radius = 50, .pos = vec2(600, 600)});

    polygon2 shape;
    shape.push_back(double2{0, 0});
    shape.push_back(double2{200, 0});
    shape.push_back(double2{200, -100});
    shape.push_back(double2{100, -100});
    shape.push_back(double2{100, -200});
    shape.push_back(double2{0, -200});
    SetShape(bodies[0], shape);

    shape.clear();
    for (int i = 0; i < 24; i++) {
	auto a = 2 * PI / 24 * i;
	shape.push_back(double2{cos(a), sin(a)} * 100);
    }
    shape.push_back(double2{0, 0});
    SetShape(bodies[1], shape);

    shape.clear();
    for (int i = 0; i < 24; i++) {
	auto a = 2 * PI / 24 * i;
	shape.push_back(double2{cos(a), sin(a)} * 50);
    }
    shape.push_back(double2{0, 0});
    SetShape(bodies[2], shape);

    mat3 ortho;
    ortho[0][0] = 2.0f / Width;
    ortho[1][1] = 2.0f / Height;
    ortho[2][2] = 1.0f;
    ortho[2][0] = -1.0f;
    ortho[2][1] = -1.0f;

    double time = 0;
    double energy = 0;
    for (Body& body : bodies) {
	energy += -glm::dot(body.pos, gravity) * body.mass + glm::dot(body.vel, body.vel) * body.mass / 2;
    }
    double energyMin = energy, energyMax = energy;

    polygon2 temp;
    vector<IContact2> contacts;

    RunEventLoop(window, [&]() {
	glClear(GL_COLOR_BUFFER_BIT);

	constexpr double dt = 0.01;
	if (gSimulate || gSimulateTick) {
	    double remaining_dt = dt;
	    while (true) {
		// resolve collisions
		// TODO Classify() call from prev iteration can tell us pairs and their contacts
		for (auto i : range(bodies.size()))
		    for (auto j : range(i + 1, bodies.size())) Interact(bodies[i], bodies[j], temp, contacts);

		SaveStates(bodies);
		Advance(bodies, remaining_dt);
		if (Classify(bodies, temp) >= 0) break;
		print("collision found!\n");

		// find collision time
		double min_dt = 0;
		double max_dt = remaining_dt;
		int i = 0;
		while (true) {
		    if (i >= 20) {
			gSimulate = false;
			print("binary search is broken (%s %s) %s\n", min_dt, max_dt, remaining_dt);
			goto exit;
		    }
		    RestoreStates(bodies);
		    double mid_dt = (min_dt + max_dt) / 2;
		    Advance(bodies, mid_dt);
		    int c = Classify(bodies, temp);
		    if (c == 0) break;
		    if (c < 0)
			max_dt = mid_dt;
		    else
			min_dt = mid_dt;
		    i += 1;
		}
		double mid_dt = (min_dt + max_dt) / 2;
		remaining_dt -= mid_dt;
	    }
	exit:

	    time += dt;
	    gSimulateTick = false;
	}

	double energy = 0;
	for (Body& body : bodies) {
	    energy += -glm::dot(body.pos, gravity) * body.mass + glm::dot(body.vel, body.vel) * body.mass / 2;
	}
	minimize(energyMin, energy);
	maximize(energyMax, energy);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	for (Body& body : bodies) {
	    vec2 p = body.pos;
	    p.x = p.x / Width * 800;
	    p.y = p.y / Height * 600;
	    monaco.render(format("pos (%f)\nvec (%f)", body.pos, body.vel), p.x, p.y, 0.3, "7FE030");
	}
	monaco.render(format("energy %s, spread %s, time %s", energy, energyMax - energyMin, time), 5, 5, 0.3,
	              "7FE030");
	glDisable(GL_BLEND);

	if (time >= 100) gSimulate = false;

	glUseProgram(shader);
	buffer.bind();
	for (const Body& body : bodies) {
	    mat3 transform(ortho);
	    transform = glm::translate(transform, vec2(body.pos));
	    transform = glm::rotate(transform, float(body.ang_pos));
	    transformUniform = transform;
	    buffer.write(cspan<vec2>(body.shapef.data(), body.shapef.size()));
	    glDrawArrays(GL_LINE_LOOP, 0, body.shapef.size());
	}
    });
    return 0;
}
