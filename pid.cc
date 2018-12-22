// auto-tunning pid controller
// - is extAcc compensation good or bad?
// - tune in presence of sensor / actuator lag
// - tune in presence of rensor and actuator signal noise
// - tune for large set of scenarios at once
// use one PID when far away AND another one when close

#include <core/format.h>
#include <core/array_deque.h>
#include <core/util.h>
#include <core/callstack.h>
#include <core/filter.h>

#include <limits>
#include <chrono>
#include <thread>
#include <random>
#include <queue>
using namespace std::chrono_literals;

#include <view/font.h>
#include <view/window.h>
#include <view/shader.h>
#include <view/vertex_buffer.h>
#include <view/glm.h>

int kSel = 0;
double kParam[3] = {-200, 0.001, -20};
bool kRecompute = false;

int gShowPos = 1;
int gShowVel = 1;
bool gShowIAcc = false;
bool gShowEAcc = true;
bool gShowAcc = false;

bool gSpring = false;
bool gPeriodic = false;
bool gConst = false;
bool gQuadDamping = false;
bool gLinearDamping = false;

void auto_tune();

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE && mods == GLFW_MOD_SHIFT)
		glfwSetWindowShouldClose(window, GL_TRUE);
	if (action == GLFW_PRESS && key == GLFW_KEY_P && mods == 0)
		gShowPos = (gShowPos + 1) % 4;
	if (action == GLFW_PRESS && key == GLFW_KEY_V && mods == 0)
		gShowVel = (gShowVel + 1) % 4;
	if (action == GLFW_PRESS && key == GLFW_KEY_I && mods == 0)
		gShowIAcc ^= 1;
	if (action == GLFW_PRESS && key == GLFW_KEY_E && mods == 0)
		gShowEAcc ^= 1;
	if (action == GLFW_PRESS && key == GLFW_KEY_A && mods == 0)
		gShowAcc ^= 1;

	if (action == GLFW_PRESS && key == GLFW_KEY_SPACE && mods == 0)
		auto_tune();

	if (action == GLFW_PRESS && key == GLFW_KEY_1 && mods == 0)
		kSel = 0;
	if (action == GLFW_PRESS && key == GLFW_KEY_2 && mods == 0)
		kSel = 1;
	if (action == GLFW_PRESS && key == GLFW_KEY_3 && mods == 0)
		kSel = 2;
	if (action == GLFW_PRESS && key == GLFW_KEY_EQUAL && mods == 0) {
		kParam[kSel] *= 1.1;
		kRecompute = true;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_MINUS && mods == 0) {
		kParam[kSel] /= 1.1;
		kRecompute = true;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_F && mods == 0) {
		kParam[kSel] *= -1;
		kRecompute = true;
	}

	if (action == GLFW_PRESS && key == GLFW_KEY_S && mods == 0) {
		gSpring ^= 1;
		kRecompute = true;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_P && mods == 0) {
		gPeriodic ^= 1;
		kRecompute = true;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_C && mods == 0) {
		gConst ^= 1;
		kRecompute = true;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_Q && mods == 0) {
		gQuadDamping ^= 1;
		kRecompute = true;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_L && mods == 0) {
		gLinearDamping ^= 1;
		kRecompute = true;
	}
}

/*
double minTimeZeroVel(double pos) {
	return 2 * sqrt(pos / maxAcc);
}

// optimal control assuming no external acc
double minTime(double pos, double vel) {
	if (pos < 0) {
		pos = -pos;
		vel = -vel;
	}

	double stopDistance = vel * vel / (2 * maxAcc);

	if (vel < 0 && stopDistance > pos) {
		// will overshoot -> brake until stop
		return -vel / maxAcc + minTimeZeroVel(stopDistance - pos);
	}

	// moving towards, but will not overshoot OR moving away
	return vel / maxAcc + minTimeZeroVel(stopDistance + pos);
}
*/

struct Params {
	double maxAcc; // m/s^2
	double minPos; // m
	double maxPosErr; // m
	double maxVelErr; // m/s
	double dt; // s
	uint inLag; // steps
	uint outLag; // steps
};

class Controller {
public:
	virtual void reset(Params params) = 0;
	virtual double execute(double pos, double vel, double time) = 0;
};

class PidController : public Controller {
public:
	double kProp = 0;
	double kInteg = 0;
	double kDeriv = 0;

	void reset(Params params) override {
		m_integral = 0;
		m_maxAcc = params.maxAcc;
	}

	double execute(double pos, double vel, double time) override {
		double prevExtAcc = (time == 0) ? 0 : ((vel - m_prevVel) / (time - m_prevTime) - m_prevAcc);
		m_prevVel = vel;
		m_prevTime = time;

		m_integral += pos;
		double acc = clamp(kProp * pos + kInteg * m_integral + kDeriv * vel /*- prevExtAcc*/, -m_maxAcc, m_maxAcc);
		m_prevAcc = acc;
		return acc;
	}

private:
	double m_integral = 0;
	double m_prevVel = 0;
	double m_prevTime = 0;
	double m_prevAcc = 0;
	double m_maxAcc = 0;
};

using ExtAcc = std::function<double(double pos, double vel, double time)>;

struct Sample {
	double pos, vel, iacc, eacc;
};

double evaluate(Controller* controller, Params params, double pos, double vel, double maxTime, bool show,
		vector<Sample>* samples, ExtAcc extAcc) {
	std::default_random_engine rnd;
	std::normal_distribution<double> normal(0, 0.0);

	LagFilter<double> posLag(params.inLag, pos);
	LagFilter<double> velLag(params.inLag, vel);
	LagFilter<double> accLag(params.outLag, 0);
	LowPassFilter<double> lowPass(1, 0);

	double stable_since = -1;

	const double dt = params.dt;
	double time = 0;
	double minPos = pos;
	uint crossings = 0;
	controller->reset(params);
	while (true) {
		if (abs(pos) <= params.maxPosErr && abs(vel) <= params.maxVelErr) {
			if (stable_since == -1) {
				stable_since = time;
			}
		} else {
			stable_since = -1;
		}

		if (time > maxTime) {
			if (show && stable_since == -1)
				print("unstable (min pos %s, crossings %s)\n", minPos, crossings);
			if (show && stable_since != -1)
				print("stabilized in %ss (min pos %s, crossings %s)\n", time, minPos, crossings);
			return (stable_since != -1) ? stable_since : std::numeric_limits<double>::infinity();
		}
		if (samples == nullptr && pos < params.minPos) {
			if (show)
				print("overshoot! (min pos %s, crossings %s)\n", minPos, crossings);
			return std::numeric_limits<double>::infinity();
		}

		double eAcc = extAcc(pos, vel, time);
		double perror = normal(rnd);
		double verror = normal(rnd);
		double aerror = normal(rnd);

		double zpos = posLag.tick(pos) + perror;
		double zvel = velLag.tick(vel) + verror;
		double zacc = controller->execute(zpos, zvel, time) + aerror;
		double iAcc = clamp(lowPass.tick(accLag.tick(zacc)), -params.maxAcc, params.maxAcc);
		double acc = eAcc + iAcc;

		if (show)
			print("time %s, pos %s, vel %s, acc %s + %s -> %s\n", time, pos, vel, eAcc, iAcc, acc);
		if (samples)
			samples->push_back({pos, vel, iAcc, eAcc});

		double p = pos;
		pos += vel * dt + acc * dt * dt / 2;
		vel += acc * dt;
		time += dt;
		minimize(minPos, pos);
		crossings += p * pos < 0;
	}
}

auto ext = [](double pos, double vel, double time) {
	double a = 0;
	if (gPeriodic)
		a += sin(time * 30) * 15;
   	if (gSpring)
		a += (10 - pos) * 10;
	if (gQuadDamping)
		a += -sign(vel) * vel * vel;
	if (gLinearDamping)
		a += -vel * 10;
	if (gConst)
		a += 60;
	return a;
};

void auto_tune() {
	std::default_random_engine rnd;
	std::normal_distribution<double> normal(0, 1);
	std::uniform_real_distribution<double> uniform(0, 1);

	Params params{.maxAcc=120, .minPos=-10, .maxPosErr=1, .maxVelErr=1, .dt=0.01, .inLag=0, .outLag=0};
	auto pid = new PidController;

	double p = 0;
	double i = 0;
	double d = 0;
	double s = std::numeric_limits<double>::infinity();
	auto start = std::chrono::system_clock::now();
	while (std::chrono::duration<double>(std::chrono::system_clock::now() - start).count() < 10) {
		double kp, ki, kd;
		kp = normal(rnd) * 10;
		ki = normal(rnd);
		kd = normal(rnd) * 2;

		pid->kProp = kp;
		pid->kInteg = ki;
		pid->kDeriv = kd;
		double ks = evaluate(pid, params, 10, 0, 100, false, nullptr, ext);

		if (ks == std::numeric_limits<double>::infinity())
			continue;
		for (int iter = 0; iter < 1000; iter++) {
			double zp, zi, zd;
			zp = kp * (1 + normal(rnd));
			zi = ki * (1 + normal(rnd));
			zd = kd * (1 + normal(rnd));

			pid->kProp = zp;
			pid->kInteg = zi;
			pid->kDeriv = zd;
			double zs = evaluate(pid, params, 10, 0, 100, false, nullptr, ext);

			if (zs < ks) {
				ks = zs;
				kp = zp;
				ki = zi;
				kd = zd;
			}
		}

		if (ks < s) {
			p = kp;
			i = ki;
			d = kd;
			s = ks;
			print("found %s\n", s);
		}
	}
	kParam[0] = p;
	kParam[1] = i;
	kParam[2] = d;
	kRecompute = true;
}

int main(int argc, char** argv) {
	InitSegvHandler();

	Params params{.maxAcc=120, .minPos=-10, .maxPosErr=0.01, .maxVelErr=0.01, .dt=0.01, .inLag=1, .outLag=1};
	auto pid = new PidController;
	pid->kProp = kParam[0];
	pid->kInteg = kParam[1];
	pid->kDeriv = kParam[2];
	vector<Sample> samples;
	evaluate(pid, params, 10, 0, 500, false, &samples, ext);

	constexpr int Width = 2550, Height = 1400;
	auto window = CreateWindow({.width=Width, .height=Height, .resizeable=false});
	glfwSetKeyCallback(window, key_callback);

	glClearColor(0.0, 0.5, 0.0, 1.0);

	Shader shader(R"END(
		#version 330 core
		uniform mat3 transform;
		uniform vec4 color;
		layout (location = 0) in vec2 pos;
	 	out vec4 vertex_color;

		void main() {
		    vec3 p = transform * vec3(pos, 1.0);
		    gl_Position = vec4(p.x, p.y, 0.0, 1.0);
			vertex_color = color;
		}

		#version 330 core
		in vec4 vertex_color;
		out vec4 color;

		void main() {
		    color = vertex_color;
		}
	)END");
	UNIFORM(mat3, transform);
	UNIFORM(vec4, color);

	VertexBuffer_vec2 buffer(500);
	std::array<vec2, 500> v;

	Font timesNewRoman("Times New Roman");
	Font arial("Arial");
	Font monaco("/System/Library/Fonts/Monaco.dfont");

	glClearColor(0, 0, 0, 1.0f);

	mat3 ortho;
	ortho[0][0] = 2.0f / Width;
	ortho[1][1] = 2.0f / Height;
	ortho[2][2] = 1.0f;
	ortho[2][0] = -1.0f;
	ortho[2][1] = -1.0f;

	RunEventLoop(window, [&]() {
		if (kRecompute) {
			samples.clear();
			pid->kProp = kParam[0];
			pid->kInteg = kParam[1];
			pid->kDeriv = kParam[2];
			evaluate(pid, params, 10, 0, 500, false, &samples, ext);
			kRecompute = false;
		}
		const uint n = min<uint>(v.size(), uint(samples.size()));

		glClear(GL_COLOR_BUFFER_BIT);

		double mx, my;
		glfwGetCursorPos(window, &mx, &my);
		my = Height - 1 - my;

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		double kTime = 5;
		double kAcc = 3;
		if (mx >= 0 && mx < Width && my >= 0 && my < Height) {
			const Sample& s = samples[round(mx / kTime)];
			auto msg = format("time %.2f, pos %.3f, vel %.2f, eacc %.2f, iacc %.2f, acc %.2f",
				mx / kTime * params.dt, s.pos, s.vel, s.eacc, s.iacc, s.eacc + s.iacc);
			timesNewRoman.render(msg, 5, 5, 0.15, "7FE030");
		}
		string_view fmt[3] = {
			"prop [%s], integ %s, deriv %s",
			"prop %s, integ [%s], deriv %s",
			"prop %s, integ %s, deriv [%s]"};
		timesNewRoman.render(format(fmt[kSel], kParam[0], kParam[1], kParam[2]), 5, 600-10, 0.15, "7FE030");

		string mods;
		if (gSpring) mods += "spring ";
		if (gPeriodic) mods += "periodic ";
		if (gQuadDamping) mods += "quad_damping ";
		if (gLinearDamping) mods += "linear_damping ";
		if (gConst) mods += "constant ";
		timesNewRoman.render(mods, 200, 600-10, 0.15, "7FE030");
		glDisable(GL_BLEND);

		glUseProgram(shader);
		buffer.bind();
		transformUniform = ortho;

		v[0] = vec2(0, Height / 2);
		v[1] = vec2(Width, Height / 2);
		buffer.write(v);
		colorUniform = vec4(1, 1, 1, 1);
		glDrawArrays(GL_LINE_STRIP, 0, 2);

		v[0] = vec2(mx, Height);
		v[1] = vec2(mx, 0);
		buffer.write(v);
		colorUniform = vec4(1, 1, 1, 1);
		glDrawArrays(GL_LINE_STRIP, 0, 2);

		if (gShowPos) {
			double kPos = 70 * (1 << (gShowPos - 1));
			for (uint i = 0; i < n; i++)
				v[i] = vec2(i * kTime, samples[i].pos * kPos + Height / 2);
			buffer.write(v);
			colorUniform = vec4(0, 0.5, 1, 1);
			glDrawArrays(GL_LINE_STRIP, 0, n);
		}

		if (gShowVel) {
			double kVel = 24 * (1 << (gShowVel - 1));
			for (uint i = 0; i < n; i++)
				v[i] = vec2(i * kTime, samples[i].vel * kVel + Height / 2);
			buffer.write(v);
			colorUniform = vec4(0, 1, 0.5, 1);
			glDrawArrays(GL_LINE_STRIP, 0, n);
		}

		if (gShowIAcc) {
			for (uint i = 0; i < n; i++)
				v[i] = vec2(i * kTime, samples[i].iacc * kAcc + Height / 2);
			buffer.write(v);
			colorUniform = vec4(1, 1, 0, 1);
			glDrawArrays(GL_LINE_STRIP, 0, n);
		}

		if (gShowEAcc) {
			for (uint i = 0; i < n; i++)
				v[i] = vec2(i * kTime, samples[i].eacc * kAcc + Height / 2);
			buffer.write(v);
			colorUniform = vec4(1, 0, 0, 1);
			glDrawArrays(GL_LINE_STRIP, 0, n);
		}

		if (gShowAcc) {
			for (uint i = 0; i < n; i++)
				v[i] = vec2(i * kTime, (samples[i].eacc + samples[i].iacc) * kAcc + Height / 2);
			buffer.write(v);
			colorUniform = vec4(1, 0.5, 0, 1);
			glDrawArrays(GL_LINE_STRIP, 0, n);
		}
	});
	return 0;
}
