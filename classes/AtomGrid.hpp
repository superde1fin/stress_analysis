#pragma once

#include <array>
#include <vector>
#include <unordered_map>
#include <map>

#include "Atom.hpp"
#include "Grid.hpp"

using namespace std;

class AtomGrid: public Grid{
    public:
        AtomGrid(array<float, 3> box, int num_splits, vector<Atom>* atoms_ptr, map<int, float> radii_mapping);
        int get_status(int key_x, int key_y, int key_z);
        void add_atom(int key_x, int key_y, int key_z, Atom* atm_ptr);
        float get_density();
    private:
        map<int, float> radii_mapping;
        float free_volume(array<int, 3> key);
        vector<Atom> get_neighbors(array<int, 3> key);
        bool all_around(array<int, 3> key);
        float median_density;
    };
