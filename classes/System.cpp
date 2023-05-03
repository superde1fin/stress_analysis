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
#include "Grid.hpp"
#include "AtomGrid.hpp"
#include "MaskGrid.hpp"

using namespace std;


System::System(array<float, 3> box, ifstream& contents, int atoms_number, array<float, 3> center, array<float, 3> box_shift){
    System::box = box;
    System::box_shift = box_shift;
    System::center = center;
    System::atoms_number = atoms_number;
    System::scan_positions(contents);
    cout << "Atoms scanned: " << atoms.size() << "\n\n";
    cout << "Modifiers detected: " << modifiers.size() << "\n\n";
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

vector<Atom> System::detect_surface(float void_volume){
    cout << "Starting surface detection\n\n";
    vector<Atom> surface_atoms;
    int prev_num_surface;
    int cur_num_surface = 0;
    int mesh_size;
    int num_splits = 2;
    array<float, 3> grid_sides;
    array<int, 3> grid_spans;
    AtomGrid* grid;
    MaskGrid* masks = new MaskGrid(System::box, num_splits);
    int debug_ctr = 0;
    float cell_volume;
    bool done = false;
    while(!done){
        prev_num_surface = cur_num_surface;
        grid = new AtomGrid(System::box, num_splits, &(System::atoms), masks, System::radii_mapping, void_volume);
        cell_volume = grid -> get_cell_volume();
        if(cell_volume >= void_volume){
            grid_sides = grid -> get_float_sides();
            grid_spans = grid -> get_int_sides();
            cout << "Sides: " << grid_sides[0] << " " << grid_sides[1] << " " << grid_sides[2] << endl;
            cout << "Spans: " << grid_spans[0] << " " << grid_spans[1] << " " << grid_spans[2] << endl;
            cout << "Average mesh density: " << grid -> get_density() << endl;
            mesh_size = grid -> get_size();;
            cout << "Atom mesh size: " << mesh_size << endl;
            surface_atoms = grid -> get_surface(masks);
            cur_num_surface = surface_atoms.size();
            num_splits *= 2;
            cout << "Prev: " << prev_num_surface << " Cur: " << cur_num_surface << endl;
            cout << "Number of masks: " << masks -> get_member_num() << endl;
            cout << "-------------------------------\n";
            }
        else{done = true;}
//        if(debug_ctr >= 1){done = true;}
        }
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
