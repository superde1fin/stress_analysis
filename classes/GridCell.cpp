#include "GridCell.hpp"

using namespace std;

//Tells if the cell is empty
bool GridCell::is_empty(){
    return GridCell::atoms.empty();
    }

//Adds an atom to an existing cell
void GridCell::add_atom(Atom* atm_ptr){
    GridCell::atoms.push_back(atm_ptr);
    }

//Returns the number of atoms in a cell
int GridCell::get_size(){
    return GridCell::atoms.size();
    }

//Cell destructor, deletes all atoms from a cell
GridCell::~GridCell(){
    GridCell::atoms.clear();
    }
