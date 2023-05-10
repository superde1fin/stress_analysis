#pragma once

#include <array>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <set>

#include "Atom.hpp"
#include "AtomGrid.hpp"
#include "MaskGrid.hpp"
#include "Molecule.hpp"

using namespace std;

class System{

    public:
        System(array<float, 3> box, ifstream& contents, int atoms_number, array<float, 3> center, array<float, 3> box_shift);

        vector<Atom> detect_surface(float void_volume);
        vector<vector<float>> calc_stresses();
        vector<vector<float>> average_stresses(float binwidth);
        float system_stress();
        void isolate_surface(array<string, 2> header, string filename);
        map<string, float> get_surface_species(int htype);
    private:
        set<int> exclude;
        Molecule* scan_molecule(Atom prev_atm, Atom atm, set<int>* exclude);
        AtomGrid* grid;
        MaskGrid* masks;
        map<int, float> radii_mapping = {{1, 1.32}, {2, 0.73}, {3, 0.53}};
        map<int, map<int, float>> cutoffs = {{1, {{1, 0}, {2, 1.8}, {3, 0}, {4, 0}}}, {2, {{2, 0}, {1, 1.8}, {3, 2.7}, {4, 1.1}}}, {3, {{3, 0}, {1, 0}, {2, 2.7}, {4, 0}}}, {4, {{4, 0}, {1, 0}, {2, 1.1}, {3, 0}}}};
        map<int, string> types = {{1, "Si"}, {2, "O"}, {3, "Na"}, {4, "H"}};
        int atoms_number;
        array<float, 3> box;
        array<float, 3> box_shift;
        array<float, 3> center;
        vector<Atom> atoms;
        vector<Atom> surface_atoms;
        vector<vector<float>> surface_stresses;
        vector<vector<float>> averaged_stresses;
        vector<Atom> modifiers;

        void scan_positions(ifstream& contents);
        vector<Atom> filter_surface(vector<Atom> unfiltered);
    };
