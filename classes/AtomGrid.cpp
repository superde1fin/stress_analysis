#include <array>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>

#include "AtomGrid.hpp"
#include "MaskGrid.hpp"
#include "Grid.hpp"
#include "GridCell.hpp"
#include "Helper.hpp"

using namespace std;

//Basic AtomGrid constructor
AtomGrid::AtomGrid(array<float, 3> box, int num_splits, vector<Atom>* atoms_ptr):Grid(box, num_splits){
    array<int, 3> key;
    array<float, 3> atom_position;
    //Put each atom into its grid cell
    for(auto it = atoms_ptr -> begin(); it != atoms_ptr -> end(); it++){
        //Calculate atom position
        atom_position = it -> get_position();
        //Calculate cell position
        key = {Helper::true_modulo((int)(atom_position[0]/AtomGrid::float_sides[0]), AtomGrid::int_sides[0]), Helper::true_modulo((int)(atom_position[1]/AtomGrid::float_sides[1]), AtomGrid::int_sides[1]), Helper::true_modulo((int)(atom_position[2]/AtomGrid::float_sides[2]), AtomGrid::int_sides[2])};
        if(!AtomGrid::grid_map.count(key)){
            GridCell* gc_ptr = new GridCell();
            AtomGrid::grid_map[key] = gc_ptr;
            }
        AtomGrid::grid_map[key] -> add_atom(&*it);
        }
    }

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
            key = {Helper::true_modulo((int)(atom_position[0]/AtomGrid::float_sides[0]), AtomGrid::int_sides[0]), Helper::true_modulo((int)(atom_position[1]/AtomGrid::float_sides[1]), AtomGrid::int_sides[1]), Helper::true_modulo((int)(atom_position[2]/AtomGrid::float_sides[2]), AtomGrid::int_sides[2])};
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
vector<Atom> AtomGrid::get_surface(MaskGrid* mask_ptr, float surface_thickness){
    vector<Atom> surface_atoms;
    vector<Atom> neighbors;
    //For each non-masked grid cell
    for(auto it = AtomGrid::grid_map.begin(); it != AtomGrid::grid_map.end(); ++it){
        //If grid cell is empty
        if((it -> second) -> is_empty()){
            //Scan closest neighbors (surface)
            neighbors = AtomGrid::get_neighbors(it -> first, surface_thickness);
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
vector<Atom> AtomGrid::get_neighbors(array<int, 3> origin_key, float surface_thickness){
    if(!surface_thickness){
        surface_thickness = 1.5*(AtomGrid::float_sides[0] + AtomGrid::float_sides[1] + AtomGrid::float_sides[2])/3;
        }
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
                        if(dist < surface_thickness){
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

vector<tuple<Atom, float>> AtomGrid::find_neighbors(Atom* atm, set<int> exclude, map<int, map<int, float>> cutoffs){
    array<float, 3> atom_position = atm -> get_position();
    int depth = 1;
    vector<tuple<Atom, float>> neighbors;
    int neigh_ctr = 0;
    int neigh_size;
    array<int, 3> key;
    float min_modifyier_dist = numeric_limits<float>::infinity();
    array<int, 3> origin_key = {(int)(atom_position[0]/AtomGrid::float_sides[0]), (int)(atom_position[1]/AtomGrid::float_sides[1]), (int)(atom_position[2]/AtomGrid::float_sides[2])};
    float max_dist = -1*numeric_limits<float>::infinity();
    float dist = 0;
    float max_cutoff = -1*numeric_limits<float>::infinity();
    for(auto it_out = cutoffs.begin(); it_out != cutoffs.end(); ++it_out){
        for(auto it_in = (it_out -> second).begin(); it_in != (it_out -> second).end(); ++it_in){
            if(it_in -> second > max_cutoff){
                max_cutoff = it_in -> second;
                }
            }
        }
    while((depth <= *max_element(AtomGrid::int_sides.begin(), AtomGrid::int_sides.begin())) && (max_dist < max_cutoff)){
        for(int i = -1*depth; i <= depth; i++){
            for(int j = -1*depth; j <= depth; j++){
                for(int k = -1*depth; k <= depth; k++){
                    key = array<int, 3>{Helper::true_modulo(origin_key[0] + i, AtomGrid::int_sides[0]), Helper::true_modulo(origin_key[1] + j, AtomGrid::int_sides[1]), Helper::true_modulo(origin_key[2] + k, AtomGrid::int_sides[2])};
                    if(AtomGrid::grid_map.count(key)){
                        for(Atom* atm_ptr : AtomGrid::grid_map[key]-> atoms){
                            dist = Helper::dist(atom_position, atm_ptr -> get_position(), AtomGrid::box);
//                            cout << atm_ptr -> get_id() << " " << dist << " " << atm_ptr -> get_type() <<  endl;
                            if(dist > max_dist){max_dist = dist;}
                            if(dist < cutoffs[atm_ptr -> get_type()][atm -> get_type()] && !Helper::element_in(atm_ptr -> get_id(), exclude) && atm -> get_id() != atm_ptr -> get_id()){
                                if((atm -> get_type() == 3 || atm -> get_type() == 4) && dist < min_modifyier_dist){
                                    min_modifyier_dist = dist;
                                    neighbors.clear();
                                    }
//                                cout << "Added\n";
                                neigh_size = neighbors.size();
                                neigh_ctr = 0;
                                while(neigh_ctr < neigh_size && get<float>(neighbors[neigh_ctr]) < dist){
                                    neigh_ctr++;
                                    }
                                neighbors.insert(next(neighbors.begin(), neigh_ctr), make_tuple(*atm_ptr, dist));
                                }
                            }
                        }
                    }
                }
            }
        depth++;
        }
    return neighbors;
    }

void AtomGrid::reset_grid(vector<Atom>* atoms_ptr, map<int, map<int, float>> cutoffs){
    array<int, 3> key;
    array<float, 3> atom_position;
    for(int i = 0; i < AtomGrid::int_sides[0]; i++){
        for(int j = 0; j < AtomGrid::int_sides[1]; j++){
            for(int k = 0; k < AtomGrid::int_sides[2]; k++){
                key = {i, j, k};
                //Create new cell
                GridCell* gc_ptr = new GridCell();
                //Add to map
                AtomGrid::grid_map[key] = gc_ptr;
                }
            }
        }
    for(auto it = atoms_ptr -> begin(); it != atoms_ptr -> end(); it++){
        //Calculate atom position
        atom_position = it -> get_position();
        if(it -> get_type() == 3 || it -> get_type() == 4){
            Atom oxygen;
            set<int> exclude;
            for(tuple<Atom, float>pair : AtomGrid::find_neighbors(&*it, exclude, cutoffs)){
                Atom oxygen = get<Atom>(pair);
                if(oxygen.get_type() == 2){
                    it -> Oid = oxygen.get_id();
                    }
                }
            }
        //Calculate cell position
        key = {Helper::true_modulo((int)(atom_position[0]/AtomGrid::float_sides[0]), AtomGrid::int_sides[0]), Helper::true_modulo((int)(atom_position[1]/AtomGrid::float_sides[1]), AtomGrid::int_sides[1]), Helper::true_modulo((int)(atom_position[2]/AtomGrid::float_sides[2]), AtomGrid::int_sides[2])};
        //Add the atom into a cell
        AtomGrid::grid_map[key] -> add_atom(&*it);
        }
    }


tuple<Atom, float> AtomGrid::find_closest(Atom* atm){
    set<int> exclude;
    return AtomGrid::find_closest(atm, exclude);
    }

//Finds closest atom in the grid
tuple<Atom, float> AtomGrid::find_closest(Atom* atm, set<int>& exclude){
    array<float, 3> atom_position = atm -> get_position();
    int depth = 1;
    tuple<Atom, float> closest;
    array<int, 3> key;
    array<int, 3> origin_key = {(int)(atom_position[0]/AtomGrid::float_sides[0]), (int)(atom_position[1]/AtomGrid::float_sides[1]), (int)(atom_position[2]/AtomGrid::float_sides[2])};
    float dist;
    float min_dist = numeric_limits<float>::infinity();
    while(depth <= *max_element(AtomGrid::int_sides.begin(), AtomGrid::int_sides.begin()) && min_dist == numeric_limits<float>::infinity()){
        for(int i = -1*depth; i <= depth; i++){
            for(int j = -1*depth; j <= depth; j++){
                for(int k = -1*depth; k <= depth; k++){
                    key = array<int, 3>{Helper::true_modulo(origin_key[0] + i, AtomGrid::int_sides[0]), Helper::true_modulo(origin_key[1] + j, AtomGrid::int_sides[1]), Helper::true_modulo(origin_key[2] + k, AtomGrid::int_sides[2])};
                    if(AtomGrid::grid_map.count(key)){
                        for(Atom* atm_ptr : AtomGrid::grid_map[key]-> atoms){
                            dist = Helper::dist(atom_position, atm_ptr -> get_position(), AtomGrid::box);
                            if(dist < min_dist && atm_ptr -> get_id() != atm -> get_id() && !Helper::element_in(atm_ptr -> get_id(), exclude)){
                                min_dist = dist;
                                closest = make_tuple(*atm_ptr,dist);
                                }
                            }
                        }
                    }
                }
            }
        depth++;
        }
    return closest;
    }
