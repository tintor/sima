#pragma once
#include <core/format.h>
#include <santorini/coord.h>

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

void Print(const Action& action) {
    std::visit(overloaded{[&](NextAction a) { print("next"); },
                          [&](PlaceAction a) { print("place:%s%s", a.dest.x(), a.dest.y()); },
                          [&](MoveAction a) { print("move:%s%s:%s%s", a.src.x(), a.src.y(), a.dest.x(), a.dest.y()); },
                          [&](BuildAction a) { print("build:%s%s:%s", a.dest.x(), a.dest.y(), a.dome ? 'D' : 'T'); }},
               action);
}

void Print(const vector<Action>& actions) {
    for (const Action& action : actions) {
        Print(action);
        print(" ");
    }
    print("\n");
}
