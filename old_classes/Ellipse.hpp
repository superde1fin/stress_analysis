#pragma once
#include <array>
#include "Atom.hpp"

using namespace std;

class Ellipse{
    public:
        Ellipse(array<float, 3> origin, float thickness, Atom atm);
        bool within(Atom atm);
        //Getters
        float get_thickness();
        float get_area();
        float get_x();
        float get_y();
        float get_z();
        float get_minor();
        float get_major();
    private:
        float x;
        float y;
        float z;
        float thickness;
        float minor;
        float major;
    };
