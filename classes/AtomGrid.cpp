#include <array>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>

#include "AtomGrid.hpp"
#include "MaskGrid.hpp"
#include "Grid.hpp"
#include "GridCell.hpp"
#include "Helper.hpp"

using namespace std;

//Atom grid constructor
AtomGrid::AtomGrid(array<float, 3> box, int num_splits, vector<Atom>* atoms_ptr, MaskGrid* mask_ptr, map<int, float> radii_mapping, float void_volume):Grid(box, num_splits){
    //Record the radii of different atoms
    AtomGrid::radii_mapping = radii_mapping;
    //Calculate grid cell volume
    AtomGrid::cell_volume = AtomGrid::float_sides[0]*AtomGrid::float_sides[1]*AtomGrid::float_sides[2];
    //Record maximum detectable void
    AtomGrid::void_volume = void_volume;
    //Check that the cell is still large enough
    if(AtomGrid::cell_volume >= void_volume){
        array<int, 3> key;
        array<float, 3> atom_position;
        array<float, 3> cell_position;
        //Loop through all possible positions of the grid
        for(int i = 0; i < AtomGrid::int_sides[0]; i++){
            for(int j = 0; j < AtomGrid::int_sides[1]; j++){
                for(int k = 0; k < AtomGrid::int_sides[2]; k++){
                    key = {i, j, k};
                    //Calculate cell position
                    cell_position = {i*AtomGrid::float_sides[0], j*AtomGrid::float_sides[1], k*AtomGrid::float_sides[2]};
                    //If cell position is not masked, add the grid cell to the grid
                    if(!(mask_ptr -> in_mask(cell_position))){
                        //Create new cell
                        GridCell* gc_ptr = new GridCell();
                        //Add to map
                        AtomGrid::grid_map[key] = gc_ptr;
                        }
                    }
                }
            }


        //Reset the mask grid
        *mask_ptr = MaskGrid(box, num_splits);
        //Put each atom into its grid cell
        for(auto it = atoms_ptr -> begin(); it != atoms_ptr -> end(); it++){
            //Calculate atom position
            atom_position = it -> get_position();
            //Calculate cell position
            key = {(int)(atom_position[0]/AtomGrid::float_sides[0]), (int)(atom_position[1]/AtomGrid::float_sides[1]), (int)(atom_position[2]/AtomGrid::float_sides[2])};
            //If the cell has been initialized, then it is not masked
            if(AtomGrid::grid_map.count(key)){
                //Add the atom into a cell
                AtomGrid::grid_map[key] -> add_atom(&*it);
                }
            else{
                //Cell not initialized, this position was masked by previous mask grid
                mask_ptr -> add_mask(key);
                }
            }
        //Calculate average number density in a cell
        AtomGrid::average_density = (atoms_ptr -> size())/AtomGrid::get_size();;
        }
    }

//Density getter
float AtomGrid::get_density(){
    return AtomGrid::average_density;
    }

//A function that detects the surface
vector<Atom> AtomGrid::get_surface(MaskGrid* mask_ptr){
    vector<Atom> surface_atoms;
    vector<Atom> neighbors;
    //For each non-masked grid cell
    for(auto it = AtomGrid::grid_map.begin(); it != AtomGrid::grid_map.end(); ++it){
        //If grid cell is empty
        if((it -> second) -> is_empty()){
            //Scan closest neighbors (surface)
            neighbors = AtomGrid::get_neighbors(it -> first);
            //Add each neighbor atom to the surface
            for(auto& neigh : neighbors){
                surface_atoms.push_back(neigh);
                }
            }
        else{
            //If cell not empty, check whether it has full neighbor cells
            if(AtomGrid::all_around(it -> first, 2)){
                //If cell has 2 full layers cells around it, it can be masked
                mask_ptr -> add_mask(it -> first);
                }
            }
        }
    //Remove duplicate atoms from the surface vector
    Helper::remove_dupl(surface_atoms);
    return surface_atoms;
    }

//Returns a vector of neighboring atoms to an empty cell
vector<Atom> AtomGrid::get_neighbors(array<int, 3> origin_key){
    float dist;
    array<int, 3> key;
    array<float, 3> center;
    vector<Atom> result;
    //Detect the center of the cell in question
    center = array<float, 3>{((float)(origin_key[0]) + (float)0.5)*(AtomGrid::float_sides[0]), ((float)origin_key[1] + (float)0.5)*(AtomGrid::float_sides[1]), ((float)origin_key[2] + (float)0.5)*(AtomGrid::float_sides[2])};
    //Cycle through cells around a given cell
    for(int i = -1; i < 2; i++){
        for(int j = -1; j < 2; j++){
            for(int k = -1; k < 2; k++){
                //Create a key for the neighboring cell
                key = array<int, 3>{Helper::true_modulo(origin_key[0] + i, AtomGrid::int_sides[0]), Helper::true_modulo(origin_key[1] + j, AtomGrid::int_sides[1]), Helper::true_modulo(origin_key[2] + k, AtomGrid::int_sides[2])};
                //If such cell exists (not masked)
                if(AtomGrid::grid_map.count(key)){
                    //For each atom in the cell
                    for(Atom* atm_ptr : AtomGrid::grid_map[key] -> atoms){
                        //Calculate the distance from an atom to the center of cell in question
                        dist = Helper::dist(atm_ptr -> get_position(), center, AtomGrid::box);
                        //If the distance is less than a specific fraction of cell dimensions, add to neighbor list
                        if(dist < 1.5*(AtomGrid::float_sides[0] + AtomGrid::float_sides[1] + AtomGrid::float_sides[2])/3){
                            result.push_back(*atm_ptr);
                            }
                        }
                    }
                }
            }
        }
    return result;
    }
//Check whether a cell is surrounded by cells full of atoms
bool AtomGrid::all_around(array<int, 3> origin_key, int depth){
    array<int, 3> key;
    array<float, 3> center;
    int neigh_ctr = 0;
    //Detect the center of the cell in question
    center = array<float, 3>{((float)(origin_key[0]) + (float)0.5)*(AtomGrid::float_sides[0]), ((float)origin_key[1] + (float)0.5)*(AtomGrid::float_sides[1]), ((float)origin_key[2] + (float)0.5)*(AtomGrid::float_sides[2])};
    //Cycle through cells around a given cell
    for(int i = -1*depth; i <= depth; i++){
        for(int j = -1*depth; j <= depth; j++){
            for(int k = -1*depth; k <= depth; k++){
                //Create a key for the neighboring cell
                key = array<int, 3>{Helper::true_modulo(origin_key[0] + i, AtomGrid::int_sides[0]), Helper::true_modulo(origin_key[1] + j, AtomGrid::int_sides[1]), Helper::true_modulo(origin_key[2] + k, AtomGrid::int_sides[2])};
                //If a cell is not masked and the number of atoms in it is more than half of the average
                if(AtomGrid::grid_map.count(key) && AtomGrid::grid_map[key] -> get_size() > 0.5*AtomGrid::average_density){
                    //Count the cell as full
                    neigh_ctr++;
                    }
                }
            }
        }
    //Returns true if all cells are counted as full
    return neigh_ctr == pow(2*depth + 1, 3);
    }

float AtomGrid::get_cell_volume(){
    return AtomGrid::cell_volume;
    }
