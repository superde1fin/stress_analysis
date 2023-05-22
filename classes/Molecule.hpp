#pragma once

#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <set>

using namespace std;

class Molecule{
    public:
        Molecule(int type, string type_str, int id);
        vector<Molecule*> get_bonds();
        Molecule* add_bond(int type, string type_str, int id, int dist);
        ~Molecule();
        bool is_atom();
        bool in_molecule(int type);
        string get_name();
        Molecule* add_bond(Molecule* ml, int dist, int rec_ctr);
        Molecule* add_bond(Molecule* ml, int dist);
        string get_type_str();
        int get_type();
        int count_atoms();
        Molecule* truncate(int depth, Molecule* parent, int charge);
        Molecule* truncate(int depth);
        tuple<int, Molecule*> closest(int type);
        tuple<int, Molecule*> closest(int type, set<int> *exclude);
        vector<int> get_ids();
        int get_id();
    private:
        map<int, int> charges = {{1, 1}, {2, -2},  {3, 1}, {4, 1}};
        Molecule* copy();
        int id;
        int charge;
        int type;
        string type_str;
        vector<Molecule*> bonds;
        vector<float> distances;
        map<int, tuple<vector<int>, int, string>> get_count(set<int>*);
    };
