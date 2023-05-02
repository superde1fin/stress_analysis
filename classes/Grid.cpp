#include <array>

#include "Grid.hpp"
#include "Helper.hpp"

using namespace std;

Grid::Grid(array<float, 3> box, int num_splits){
    Grid::box = box;
    tuple<array<float, 3>, array<int, 3>> sides_tuple = Helper::calc_cell_spans(box, num_splits);
    Grid::float_sides = get<array<float, 3>>(sides_tuple);
    Grid::int_sides = get<array<int, 3>>(sides_tuple);
    }

array<int, 3> Grid::get_int_sides(){
    return Grid::int_sides;
    }

array<float, 3> Grid::get_float_sides(){
    return Grid::float_sides;
    }

Grid::~Grid(){
    for(auto it : Grid::grid_map){
        delete it.second;
        }
    }

int Grid::get_size(){
    return Grid::int_sides[0]*Grid::int_sides[1]*Grid::int_sides[2];
    }

int Grid::get_member_num(){
    return Grid::grid_map.size();
    }
