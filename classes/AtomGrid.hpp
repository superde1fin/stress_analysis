#pragma once

#include <array>
#include <vector>
#include <unordered_map>
#include <map>
#include <tuple>
#include <set>

#include "Atom.hpp"
#include "Grid.hpp"
#include "MaskGrid.hpp"

using namespace std;

class AtomGrid: public Grid{
    public:
        AtomGrid(array<float, 3> box, int num_splits, vector<Atom>* atoms_ptr, MaskGrid* mask_ptr, map<int, float> radii_mapping, float void_volume);
        AtomGrid(array<float, 3> box, int num_splits, vector<Atom>* atoms_ptr);
        tuple<Atom, float> find_closest(Atom* atm);
        tuple<Atom, float> find_closest(Atom* atm, set<int>& exclude);
        int get_status(int key_x, int key_y, int key_z);
        void add_atom(int key_x, int key_y, int key_z, Atom* atm_ptr);
        float get_density();
        vector<Atom> get_surface(MaskGrid* mask_ptr, float surface_thickness);
        float get_cell_volume();
        void reset_grid(vector<Atom>* atoms_ptr, map<int, map<int, float>> cutoffs);
    private:
        vector<tuple<Atom, float>> find_neighbors(Atom* atm, set<int> exclude, map<int, map<int, float>> cutoffs);
        float cell_volume;
        float void_volume;
        float average_density;
        map<int, float> radii_mapping;
        float free_volume(array<int, 3> key);
        vector<Atom> get_neighbors(array<int, 3> key, float surface_thickness);
        bool all_around(array<int, 3> key, int depth);
        float median_density;
    friend class System;
    };
