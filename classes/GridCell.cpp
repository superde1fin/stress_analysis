#include "GridCelll.hpp"

using namespace std;

bool GridCelll::is_empty(){
    return GridCelll::atoms.empty();
    }

void GridCelll::add_atom(Atom* atm_ptr){
    GridCelll::atoms.push_back(atm_ptr);
    }

int GridCelll::get_size(){
    return GridCelll::atoms.size();
    }

GridCelll::~GridCelll(){
    GridCelll::atoms.clear();
    }
