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

AtomGrid::AtomGrid(array<float, 3> box, int num_splits, vector<Atom>* atoms_ptr, MaskGrid* mask_ptr, map<int, float> radii_mapping):Grid(box, num_splits){
    AtomGrid::radii_mapping = radii_mapping;
    array<int, 3> key;
    array<float, 3> atom_position;
    array<float, 3> cell_position;
    for(int i = 0; i < AtomGrid::int_sides[0]; i++){
        for(int j = 0; j < AtomGrid::int_sides[1]; j++){
            for(int k = 0; k < AtomGrid::int_sides[2]; k++){
                key = {i, j, k};
                cell_position = {i*AtomGrid::float_sides[0], j*AtomGrid::float_sides[1], k*AtomGrid::float_sides[2]};
                if(!(mask_ptr -> in_mask(cell_position))){
                    GridCell* gc_ptr = new GridCell();
                    AtomGrid::grid_map[key] = gc_ptr;
                    }
                }
            }
        }


    *mask_ptr = MaskGrid(box, num_splits);
    for(auto it = atoms_ptr -> begin(); it != atoms_ptr -> end(); it++){
        atom_position = it -> get_position();
        key = {(int)(atom_position[0]/AtomGrid::float_sides[0]), (int)(atom_position[1]/AtomGrid::float_sides[1]), (int)(atom_position[2]/AtomGrid::float_sides[2])};
        if(AtomGrid::grid_map.count(key)){
            AtomGrid::grid_map[key] -> add_atom(&*it);
            }
        else{
            mask_ptr -> add_mask(key);
            }
        }
    AtomGrid::average_density = (atoms_ptr -> size())/AtomGrid::get_size();;
    }


float AtomGrid::get_density(){
    return AtomGrid::average_density;
    }


vector<Atom> AtomGrid::get_surface(MaskGrid* mask_ptr, float void_volume){
    vector<Atom> surface_atoms;
    vector<Atom> neighbors;
    array<int, 3>test_key = {15, 126, 88};
    cout << AtomGrid::grid_map.count(test_key) << endl;
    for(auto it = AtomGrid::grid_map.begin(); it != AtomGrid::grid_map.end(); ++it){
//        for(auto* neigh : (it -> second) -> atoms){
//            surface_atoms.push_back(*neigh);
//            }
        
        cout << "Key: " << (it -> first)[0] << " " << (it -> first)[1] << " " << (it -> first)[2] << endl;
//        cout << "Cell empty: " << (it -> second) -> is_empty() << endl;
        if((it -> second) -> is_empty()){
            neighbors = AtomGrid::get_neighbors(it -> first);
            for(auto& neigh : neighbors){
                surface_atoms.push_back(neigh);
                }
            }
        else{
            if(AtomGrid::all_around(it -> first, void_volume)){
                mask_ptr -> add_mask(it -> first);
                }
            }
        }
//    Helper::remove_dupl(surface_atoms);
    return surface_atoms;
    }


vector<Atom> AtomGrid::get_neighbors(array<int, 3> origin_key){
    float dist;
    array<int, 3> key;
    array<float, 3> center;
    vector<Atom> result;
    for(int i = -1; i < 2; i++){
        for(int j = -1; j < 2; j++){
            for(int k = -1; k < 2; k++){
                key = array<int, 3>{Helper::true_modulo(origin_key[0] + i, AtomGrid::int_sides[0]), Helper::true_modulo(origin_key[1] + j, AtomGrid::int_sides[1]), Helper::true_modulo(origin_key[2] + k, AtomGrid::int_sides[2])};
                center = array<float, 3>{((float)(origin_key[0]) + (float)0.5)*(AtomGrid::float_sides[0]), ((float)origin_key[1] + (float)0.5)*(AtomGrid::float_sides[1]), ((float)origin_key[2] + (float)0.5)*(AtomGrid::float_sides[2])};
                if(AtomGrid::grid_map[key]){
                    for(Atom* atm_ptr : AtomGrid::grid_map[key] -> atoms){
                        dist = Helper::dist(atm_ptr -> get_position(), center, AtomGrid::box);
                        if(dist < (AtomGrid::float_sides[0] + AtomGrid::float_sides[1] + AtomGrid::float_sides[2])/3){
                            result.push_back(*atm_ptr);
                            }
                        }
                    }
                }
            }
        }
    return result;
    }

bool AtomGrid::all_around(array<int, 3> origin_key, float void_volume){
    array<int, 3> key;
    array<float, 3> center;
    int neigh_ctr = 0;
    for(int i = -1; i < 2; i++){
        for(int j = -1; j < 2; j++){
            for(int k = -1; k < 2; k++){
                key = array<int, 3>{Helper::true_modulo(origin_key[0] + i, AtomGrid::int_sides[0]), Helper::true_modulo(origin_key[1] + j, AtomGrid::int_sides[1]), Helper::true_modulo(origin_key[2] + k, AtomGrid::int_sides[2])};
                center = array<float, 3>{((float)(origin_key[0]) + (float)0.5)*(AtomGrid::float_sides[0]), ((float)origin_key[1] + (float)0.5)*(AtomGrid::float_sides[1]), ((float)origin_key[2] + (float)0.5)*(AtomGrid::float_sides[2])};
                if(AtomGrid::grid_map[key] && AtomGrid::free_volume(key) < void_volume){
                    neigh_ctr++;
                    }
                }
            }
        }
    return neigh_ctr == 27;
    }

float AtomGrid::free_volume(array<int, 3> key){
    float cell_volume = AtomGrid::float_sides[0]*AtomGrid::float_sides[1]*AtomGrid::float_sides[2];
    float atom_volume = 0;
    for(Atom* atm_ptr : AtomGrid::grid_map[key] -> atoms){
        atom_volume += (4/3)*M_PI*pow(AtomGrid::radii_mapping[atm_ptr -> get_type()], 3);
        }
    return cell_volume - atom_volume;
    }
