#pragma once

#include <array>
#include <vector>

#include "Atom.hpp"

class Bubble{
    public:
        Bubble(array<float, 3> position, array<float, 3> span);
        bool is_void(vector<Atom> system_atoms, array<float, 3> box);
        Atom get_closest();
        array<int, 3> get_mapping();
        array<float, 3> get_position();
        array<float, 3> get_span();
    private:
        Atom closest;
        float x;
        float y;
        float z;
        float x_span;
        float y_span;
        float z_span;
        array<float, 3> position;
        array<float, 3> span;
    };
