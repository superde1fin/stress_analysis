#pragma once

#include <array>
#include <fstream>
#include <vector>
#include <map>
#include <limits>

#include "Atom.hpp"
#include "Ellipse.hpp"
#include "Sector.hpp"

using namespace std;

class System{

    public:
        System(array<float, 3> box, ifstream& contents, int atoms_number, array<float, 3> center, array<float, 3> box_shift);

        vector<Atom> detect_surface();
        vector<vector<float>> calc_stresses();
        vector<vector<float>> average_stresses(float binwidth);
        float system_stress();
        void isolate_surface(array<string, 2> header, string filename);
    private:
        map<int, float> radii_mapping = {{1, 1}, {2, 1}, {3, 1}};
        int atoms_number;
        array<float, 3> box;
        array<float, 3> box_shift;
        array<float, 3> center;
        vector<Atom> atoms;
        vector<Atom> surface_atoms;
        vector<vector<float>> surface_stresses;
        vector<vector<float>> averaged_stresses;
        vector<Atom> modifiers;
        float min_radius = numeric_limits<float>::infinity();

        void scan_positions(ifstream& contents);
        vector<vector<Atom>> split_atoms(float step);
};
