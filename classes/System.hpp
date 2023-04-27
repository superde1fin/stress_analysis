#pragma once

#include <array>
#include <fstream>
#include <vector>
#include <map>
#include <limits>

#include "Atom.hpp"

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
        float mesh_density;
        vector<array<array<float, 2>, 3>> grid_masks;
        array<float, 3> grid_sizes;
        map<int, float> radii_mapping = {{1, 1.32}, {2, 0.73}, {3, 0.53}};
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
        vector<vector<vector<vector<Atom>>>> split_atoms();
        vector<Atom> iter_get_surface(vector<vector<vector<vector<Atom>>>>& mesh);
        vector<Atom> mesh_closest(vector<vector<vector<vector<Atom>>>>& mesh, array<int, 3> origin_key, bool& all_around);
        vector<Atom> filter_surface(vector<Atom> unfiltered);
};
