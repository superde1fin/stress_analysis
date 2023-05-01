#include <array>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <limits>
#include <map>
#include <chrono>

#include "System.hpp"
#include "Atom.hpp"
#include "Helper.hpp"

using namespace std;
using namespace std::chrono;


System::System(array<float, 3> box, ifstream& contents, int atoms_number, array<float, 3> center, array<float, 3> box_shift){
    System::box = box;
    System::box_shift = box_shift;
    System::center = center;
    System::atoms_number = atoms_number;
    System::scan_positions(contents);
    cout << "Atoms scanned: " << atoms.size() << "\n\n";
    cout << "Modifiers detected: " << modifiers.size() << "\n\n";
    }


vector<vector<vector<vector<Atom>>>> System::split_atoms(){
    vector<vector<vector<vector<Atom>>>> mesh;
    int x_span = round(System::box[0]/System::grid_sizes[0]);
    int y_span = round(System::box[1]/System::grid_sizes[1]);
    int z_span = round(System::box[2]/System::grid_sizes[2]);
    array<float, 3> position;
    array<float, 3> atom_position;
    for(int i = 0; i < x_span; i++){
        vector<vector<vector<Atom>>> y_vect(y_span);
        mesh.push_back(y_vect);
        for(int j = 0; j < y_span; j++){
            vector<vector<Atom>> z_vect(z_span);
            mesh[i][j] = z_vect;
            for(int k = 0; k < z_span; k++){
                vector<Atom> atom_vect;
                position = array<float, 3>{i*System::grid_sizes[0], j*System::grid_sizes[1], k*System::grid_sizes[2]};
                mesh[i][j][k]= atom_vect;
                }
            }
        }

    for(Atom& atm : System::atoms){
        atom_position = atm.get_position();
        mesh[(int)(atom_position[0]/System::grid_sizes[0])][(int)(atom_position[1]/System::grid_sizes[1])][(int)(atom_position[2]/System::grid_sizes[2])].push_back(atm);
        }
    vector<float> densities;
    for(int i = 0; i < x_span; i++){
        for(int j = 0; j < y_span; j++){
            for(int k = 0; k < z_span; k++){
                densities.push_back(mesh[i][j][k].size());
                }
            }
        }
    sort(densities.begin(), densities.end());
    int dens_size = densities.size();
    System::mesh_density = dens_size%2 ? (densities[(int)(dens_size/2)] + densities[(int)(dens_size/2) + 1])/2 : densities[(int)(dens_size/2)];
    return mesh;
    }

vector<Atom> System::iter_get_surface(vector<vector<vector<vector<Atom>>>>& mesh){
    vector<Atom> surface_atoms;
    vector<Atom> closest;
    bool all_around;
    bool ever_empty = false, in_mask;
    array<int, 3> key;
    int i, j, k;
    int x_span = round(System::box[0]/System::grid_sizes[0]);
    int y_span = round(System::box[1]/System::grid_sizes[1]);
    int z_span = round(System::box[2]/System::grid_sizes[2]);
    Atom atm;
    array<float, 3> position;
    for(i = 0; i < x_span; i++){
        for(j = 0; j < y_span; j++){
            for(k = 0; k < z_span; k++){
                in_mask = false;
                position = array<float, 3>{i*System::grid_sizes[0], j*System::grid_sizes[1], k*System::grid_sizes[2]};
                for(auto& mask : System::grid_masks){
                    if(mask[0][0] <= position[0] && position[0] < mask[0][1] && mask[1][0] <= position[1] && position[1] < mask[1][1] && mask[2][0] <= position[2] && position[2] < mask[2][1]){
                        in_mask = true;
                        break;
                        }
                    }
                if(!in_mask){
                    key = array<int, 3>{i, j, k};
                    all_around = false;
                    closest = System::mesh_closest(mesh, key, all_around);
                    if(all_around){
                        System::grid_masks.push_back(array<array<float, 2>, 3>{array<float, 2>{i*System::grid_sizes[0], (i + 1)*System::grid_sizes[0]}, array<float, 2>{j*System::grid_sizes[1], (j + 1)*System::grid_sizes[1]}, array<float, 2>{k*System::grid_sizes[2], (k + 1)*System::grid_sizes[2]}});
                        }
                    else{
                        ever_empty = true;
                        for(Atom& atm : closest){
//                        for(Atom& atm : mesh[i][j][k]){
                            surface_atoms.push_back(atm);
                            }
                        }
                    }
                }
            }
        }
    if(!ever_empty){
        System::grid_masks.clear();
        }
    Helper::remove_dupl(surface_atoms);
    return surface_atoms;
    }

vector<Atom> System::mesh_closest(vector<vector<vector<vector<Atom>>>>& mesh, array<int, 3> origin_key, bool& all_around){
    float dist;
    int i, j, k;
    array<int, 3> key;
    array<float, 3> center;
    vector<Atom> result;
    int x_span = round(System::box[0]/System::grid_sizes[0]);
    int y_span = round(System::box[1]/System::grid_sizes[1]);
    int z_span = round(System::box[2]/System::grid_sizes[2]);
    int neigh_ctr = 0;
    for(i = -1; i < 2; i++){
        for(j = -1; j < 2; j++){
            for(k = -1; k < 2; k++){
                key = array<int, 3>{Helper::true_modulo(origin_key[0] + i, x_span), Helper::true_modulo(origin_key[1] + j, y_span), Helper::true_modulo(origin_key[2] + k, z_span)};
                center = array<float, 3>{((float)(origin_key[0]) + (float)0.5)*(System::grid_sizes[0]), ((float)origin_key[1] + (float)0.5)*(System::grid_sizes[1]), ((float)origin_key[2] + (float)0.5)*(System::grid_sizes[1])};
                for(Atom& atm : mesh[key[0]][key[1]][key[2]]){
                    dist = Helper::dist(atm.get_position(), center, System::box);
//                    if(dist < (System::grid_sizes[0] + System::grid_sizes[1] + System::grid_sizes[2])/6){
                    if(1){
                        result.push_back(atm);
                        }
                    }
                if(mesh[key[0]][key[1]][key[2]].size() > 0.3*System::mesh_density){
//                if(!mesh[key[0]][key[1]][key[2]].empty()){
                    neigh_ctr++;
                    }
                }
            }
        }
    if(neigh_ctr == 27){
        all_around = true;
        }
    return result;
    }

vector<Atom> System::filter_surface(vector<Atom> unfiltered){
    vector<Atom> filtered;
    tuple<Atom, float> closest;
    for(Atom& atm : unfiltered){
        closest = Helper::find_closest(atm, unfiltered, System::box);
        if(get<Atom>(closest).get_type() != atm.get_type() && get<float>(closest) <= 1.3*(System::radii_mapping[atm.get_type()] + System::radii_mapping[get<Atom>(closest).get_type()])){
            filtered.push_back(atm);
            }
        }
    return filtered;
    }

vector<Atom> System::detect_surface(){
    cout << "Starting surface detection\n\n";
    vector<Atom> surface_atoms;
    int prev_num_surface;
    int cur_num_surface = 0;
    int mesh_size;
    int num_splits = 2;
    vector<vector<vector<vector<Atom>>>> mesh;
    do{
        prev_num_surface = cur_num_surface;
        System::grid_sizes = Helper::calc_cell_spans(System::box, num_splits);
        cout << "Spans: " << System::grid_sizes[0] << " " << System::grid_sizes[1] << " " << System::grid_sizes[2] << endl;
        mesh = System::split_atoms();
        cout << "Average mesh density: " << System::mesh_density << endl;
        mesh_size = mesh.size()*mesh[0].size()*mesh[0][0].size();
        cout << "Atom mesh size: " << mesh_size << endl;
        surface_atoms = iter_get_surface(mesh);
        //surface_atoms = System::filter_surface(surface_atoms);
        cur_num_surface = surface_atoms.size();
        num_splits *= 2;
        cout << "Prev: " << prev_num_surface << " Cur: " << cur_num_surface << endl;
        cout << "-------------------------------\n";
        }while(!System::grid_masks.size());
//        }while(mesh_size < System::atoms_number || cur_num_surface == 0);
//        }while(cur_num_surface > prev_num_surface || cur_num_surface == 0);
//        }while((abs(prev_num_surface - cur_num_surface) >= 0.05*cur_num_surface || cur_num_surface == 0) && mesh_size < System::atoms_number);
    System::surface_atoms = surface_atoms;
    return System::surface_atoms;
    }

void System::scan_positions(ifstream& contents){
    string line;
    vector<string> split_line;
    for(int i = 0; i < System::atoms_number; i++){
        getline(contents, line);
        split_line = Helper::split(line, " ");
        Atom* atm = new Atom();
        atm -> set_pure_line(line);
        atm -> set_id(stoi(split_line[0]));
        atm -> set_type(stoi(split_line[1]));
        atm -> set_position(array<float, 3>{fmod(stof(split_line[2]) - System::box_shift[0], System::box[0]), fmod(stof(split_line[3]) - System::box_shift[1], System::box[1]), fmod(stof(split_line[4]) - System::box_shift[2], System::box[2])});
        atm -> set_stress_tensor(array<float, 6>{stof(split_line[5]), stof(split_line[6]), stof(split_line[7]),stof(split_line[8]), stof(split_line[9]), stof(split_line[10])});
        atm -> set_radius(System::radii_mapping[atm -> get_type()]);
        System::atoms.push_back(*atm);
//        cout << atm -> get_z() << " " << System::box_shift[2] << endl << line << endl;
        if(atm ->get_radius() < System::min_radius){System::min_radius = atm -> get_radius();}
        if(atm -> get_type() != 1 && atm -> get_type() != 2 ){
            System::modifiers.push_back(*atm);
            }
        delete atm;
        }

    }

vector<vector<float>> System::calc_stresses(){
    cout << "Beginning the surface stress calculations\n\n";
    tuple<Atom, float> closest;
    vector<vector<float>> stress_function;
    vector<float> tuple;
    if(System::surface_atoms.empty()){
        System::detect_surface();
        }
    for(Atom& atm : System::surface_atoms){
        if(atm.get_type() == 1){
            closest = Helper::find_closest(atm, System::modifiers, System::box);
            tuple.clear();
            tuple.push_back(get<float>(closest));
            tuple.push_back(atm.get_ave_stress());
            stress_function.push_back(tuple);
            }
        }
    System::surface_stresses = stress_function;
    return stress_function;
    }

vector<vector<float>> System::average_stresses(float binwidth){
    sort(System::surface_stresses.begin(), System::surface_stresses.end(), [=](vector<float>& vect1, vector<float>& vect2){return vect1[0] < vect2[0];});
    float temp_stress_value = 0.0;
    int datapoint_number = 0, group_number = (int)(System::surface_stresses[0][0]/binwidth);
    float prev_dist = -1.0;

    vector<vector<float>> averaged_stresses;
    vector<float> tuple;

    int stresses_vect_size = System::surface_stresses.size();
    for(int i = 0; i < stresses_vect_size; i++){
        if(prev_dist < group_number*binwidth && System::surface_stresses[i][0] >= group_number*binwidth){
            tuple.clear();
            tuple.push_back(group_number*binwidth);
            tuple.push_back(datapoint_number == 0 ? 0 : temp_stress_value/datapoint_number);
            averaged_stresses.push_back(tuple);
            group_number++;
            temp_stress_value = 0;
            datapoint_number = 0;
            }
        temp_stress_value += System::surface_stresses[i][1];
        datapoint_number++;
        prev_dist = System::surface_stresses[i][0];
        }
        tuple.clear();
        tuple.push_back(group_number*binwidth);
        tuple.push_back(datapoint_number == 0 ? 0 : temp_stress_value/datapoint_number);
        averaged_stresses.push_back(tuple);

    System::averaged_stresses = averaged_stresses;
    return averaged_stresses;
    }

float System::system_stress(){
    float total_stress = 0;
    for(Atom& atm : System::atoms){
        total_stress += atm.get_ave_stress();
        }
    return total_stress/System::atoms_number;
    }

void System::isolate_surface(array<string, 2> header, string filename){
    string atoms_num = Helper::to_str((int)System::surface_atoms.size());
    string contents = header[0] + atoms_num + "\n" + header[1];
    ofstream myfile;
    myfile.open(filename);
    myfile << contents;
    for(Atom& atm : System::surface_atoms){
        myfile << (atm.get_pure_line() + "\n");
        }
    myfile.close();
    }
