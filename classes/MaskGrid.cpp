#include <array>

#include "MaskGrid.hpp"
#include "Grid.hpp"
#include "GridCell.hpp"

using namespace std;

MaskGrid::MaskGrid(array<float, 3> box, int num_splits):Grid(box, num_splits){}

//Adds a masked grid cell
void MaskGrid::add_mask(array<int, 3> key){
    //If not masked already
    if(!MaskGrid::grid_map[key]){
        //Create new cell
        GridCell* gc_ptr = new GridCell();
        //Add to maksk
        MaskGrid::grid_map[key] = gc_ptr;
        }
    }

//Checks whether a specific key is in mask
bool MaskGrid::in_mask(array<float, 3> position){
    return MaskGrid::grid_map[array<int, 3>{(int)(position[0]/MaskGrid::float_sides[0]), (int)(position[1]/MaskGrid::float_sides[1]), (int)(position[2]/MaskGrid::float_sides[2])}];
    }
