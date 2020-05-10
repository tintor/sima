#pragma once

#include <core/std.h>
#include <geom/convex_body.h>
#include <geom/quaternion.h>

// shapes can be concave composites, but composed of convex pieces
struct ConvexShape {
    cmesh3 mesh;

    double elasticity;
    double drag;  // TODO should depend on orientation
    double staticFriction;
    double dynamicFriction;
};

struct State {
    double3 position;
    quat orientation;
    double3 velocity;
    double3 rotation;
};

struct Body {
    vector<ConvexShape> shapes;
    double radius;

    double imass;   // inverted mass (=0 if non-interactive)
    double33 imoi;  // inverted moment of intertia (=0 if non-interactive)

    State present;
    State future;

    double33 orientationMat;
};

struct Joint {
    Body* bodyA;
    Body* bodyB;

    // anchors are in body frame
    double3 anchorA;
    double3 anchorB;

    // extra for hinge joint
    bool isHinge = false;
    // axises are in body frame
    double3 axisA;
    double3 axisB;
    // TODO range of motion: angle restrictions
    // TODO motor: applying force through hinge
};

struct System {
    double time;
    vector<unique_ptr<Body>> bodies;
    vector<unique_ptr<Joint>> joints;
};

// reads Body::present, and writes to Body::future
void KinematicStep(System& system, double dtime);

void Apply(System& system);

struct Contact {
    Body* bodyA;
    Body* bodyB;
    ConvexShape* shapeA;
    ConvexShape* shapeB;
    double3 normal;
    vector<double3> points;
};

// returns true if changes were made
bool Resolve(const Contact& contact);

void ResolveConstraints(System& system);

// finds time of first collision from present state to future state during dtime
double FindCollisionTime(System& system, double dtime);

void Advance(System& system, double dtime);
