#pragma once

#include <vector>

#include "Atom.hpp"

using namespace std;

class GridCelll{
    public:
        GridCelll();
        ~GridCelll();
        bool is_empty();
        int get_size();
    private:
        vector<Atom*> atoms;
        void add_atom(Atom* atm_ptr);
    friend class Grid;
    };
