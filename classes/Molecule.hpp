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
        Molecule* add_bond(int type, string type_str, int id);
        ~Molecule();
        bool is_atom();
        bool in_molecule(int type);
        string get_name();
        Molecule* add_bond(Molecule* ml);
        string get_type_str();
        int get_type();
        int count_atoms();
        Molecule* truncate(int depth, Molecule* parent);
        tuple<int, Molecule*> closest(int type);
        tuple<int, Molecule*> closest(int type, set<int> *exclude);
        vector<int> get_ids();
    private:
        Molecule* copy();
        int id;
        int charge;
        int type;
        string type_str;
        vector<Molecule*> bonds;
        map<int, tuple<vector<int>, int, string>> get_count(set<int>*);
    };
