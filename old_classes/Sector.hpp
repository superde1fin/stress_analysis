#pragma once
#include <vector>

#include "Ellipse.hpp"
#include "Atom.hpp"

class Sector{
    public:
        Sector(Ellipse crcl, Atom atm, bool* success_indicator);
        bool in_segments(Ellipse crcl, Atom atm);
        static float get_area(Ellipse crcl);
        static void reset();

    private:
        static vector<array<float, 2>> segment_union;
        float lower_bound;
        float upper_bound;
        float thickness;

        float get_half_angle(Ellipse crcl, Atom atm);
        float get_angular_position(Ellipse crcl, Atom atm);
        void add2union();
    };

