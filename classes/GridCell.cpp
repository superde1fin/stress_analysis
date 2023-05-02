#include "GridCell.hpp"

using namespace std;

bool GridCell::is_empty(){
    return GridCell::atoms.empty();
    }

void GridCell::add_atom(Atom* atm_ptr){
    GridCell::atoms.push_back(atm_ptr);
    }

int GridCell::get_size(){
    return GridCell::atoms.size();
    }

GridCell::~GridCell(){
    GridCell::atoms.clear();
    }
