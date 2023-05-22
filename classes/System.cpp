#include <array>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <limits>
#include <map>
#include <set>

#include "System.hpp"
#include "Atom.hpp"
#include "Helper.hpp"
#include "Grid.hpp"
#include "AtomGrid.hpp"
#include "MaskGrid.hpp"

using namespace std;


//System class constructor
System::System(array<float, 3> box, ifstream& contents, int atoms_number, array<float, 3> center, array<float, 3> box_shift){
    System::box = box;
    System::box_shift = box_shift;
    System::center = center;
    System::atoms_number = atoms_number;
    //Scans all atom positions into atom objects
    System::scan_positions(contents);
    cout << "Atoms scanned: " << System::atoms.size() << "\n\n";
    cout << "Modifiers detected: " << System::modifiers.size() << "\n\n";
    }

//Regulates the surface detection performed by the AtomGrid class
vector<Atom> System::detect_surface(float void_volume){
    cout << "Starting surface detection\n\n";
    vector<Atom> surface_atoms;
    int cur_num_surface = 0;
    int mesh_size;
    int num_splits = 2;
    array<float, 3> grid_sides;
    array<int, 3> grid_spans;
    //Create an empty mask grid
    System::masks = new MaskGrid(System::box, num_splits);
    float cell_volume;
    bool done = false;
    AtomGrid* prev_grid;
    //Loop while the grid volume is greater than minimum detected void
    while(!done){
        prev_grid = System::grid;
        //Create the AtomGrid object
        System::grid = new AtomGrid(System::box, num_splits, &(System::atoms), System::masks, System::radii_mapping, void_volume);
        //Calculate the grid cell olume
        cell_volume = System::grid -> get_cell_volume();
        //If volume is large enough, continue
        if(cell_volume >= void_volume){
            //Data for console output
            grid_sides = System::grid -> get_float_sides();
            grid_spans = System::grid -> get_int_sides();
            cout << "Sides: " << grid_sides[0] << " " << grid_sides[1] << " " << grid_sides[2] << endl;
            cout << "Spans: " << grid_spans[0] << " " << grid_spans[1] << " " << grid_spans[2] << endl;
            cout << "Average mesh density: " << grid -> get_density() << endl;
            mesh_size = System::grid -> get_size();
            cout << "Atom mesh size: " << mesh_size << endl;
            //Run the grid surface detection function
            surface_atoms = System::grid -> get_surface(System::masks);
            cur_num_surface = surface_atoms.size();
            //Increase the rastering definition
            num_splits *= 2;
            cout << "Current number of surface atoms: " << cur_num_surface << endl;
            cout << "Number of masks: " << masks -> get_member_num() << endl;
            cout << "-------------------------------\n";
            }
        else{done = true;}
        }
    System::grid = prev_grid;
    //Record the acquired surface as a system member
    System::surface_atoms = surface_atoms;
    return System::surface_atoms;
    }

//A function that parses the initial .dump file into Atom objects
void System::scan_positions(ifstream& contents){
    string line;
    vector<string> split_line;
    int i;
    //Cycle through lines
    for(i = 0; i < System::atoms_number && getline(contents, line); i++){
        //Split the line by spaces
        split_line = Helper::split(line, " ");
        //Create new Atom object
        Atom* atm = new Atom();
        //Parse the position and stresses
        atm -> set_pure_line(line);
        atm -> set_id(stoi(split_line[0]));
        atm -> set_type(stoi(split_line[1]));
        atm -> set_position(array<float, 3>{fmod(stof(split_line[2]) - System::box_shift[0], System::box[0]), fmod(stof(split_line[3]) - System::box_shift[1], System::box[1]), fmod(stof(split_line[4]) - System::box_shift[2], System::box[2])});
        atm -> set_stress_tensor(array<float, 6>{stof(split_line[5]), stof(split_line[6]), stof(split_line[7]),stof(split_line[8]), stof(split_line[9]), stof(split_line[10])});
        atm -> set_radius(System::radii_mapping[atm -> get_type()]);
        //Add the atom to the system
        System::atoms.push_back(*atm);
        //Record the modifier atoms in the modifiers vector
        if(atm -> get_type() != 1 && atm -> get_type() != 2 ){
            System::modifiers.push_back(*atm);
            }
        delete atm;
        }
    //A safe guard from incorrect atom number listed in .dump
    if(i < System::atoms_number){
        System::atoms_number = i;
        }
    }

//A function that calculates per atom average stress as a function to the closest modifier
vector<vector<float>> System::calc_stresses(){
    cout << "Beginning the surface stress calculations\n\n";
    tuple<Atom, float> closest;
    vector<vector<float>> stress_function;
    vector<float> tuple;
    //For each modifier atom
    for(Atom& atm : System::surface_atoms){
        //Checking that modifiers recorded correctly
        if(atm.get_type() == 1){
            //Get closest atom
            closest = Helper::find_closest(atm, System::modifiers, System::box);
            tuple.clear();
            //Put the distance and stress into a resulting vector
            tuple.push_back(get<float>(closest));
            tuple.push_back(atm.get_ave_stress());
            stress_function.push_back(tuple);
            }
        }
    //Record the stresses vector as a system member
    System::surface_stresses = stress_function;
    return stress_function;
    }

//Average the stresses by distance regions
vector<vector<float>> System::average_stresses(float binwidth){
    //Sort the stresses vector by distance to closest modifier
    sort(System::surface_stresses.begin(), System::surface_stresses.end(), [=](vector<float>& vect1, vector<float>& vect2){return vect1[0] < vect2[0];});
    float temp_stress_value = 0.0;
    int datapoint_number = 0, group_number = (int)(System::surface_stresses[0][0]/binwidth);
    float prev_dist = -1.0;

    vector<vector<float>> averaged_stresses;
    vector<float> tuple;

    int stresses_vect_size = System::surface_stresses.size();
    //For each dist-stress pair
    for(int i = 0; i < stresses_vect_size; i++){
        //If in the next distance region
        if(prev_dist < group_number*binwidth && System::surface_stresses[i][0] >= group_number*binwidth){
            tuple.clear();
            //Record average distance
            tuple.push_back(group_number*binwidth);
            //Record average stress
            tuple.push_back(datapoint_number == 0 ? 0 : temp_stress_value/datapoint_number);
            averaged_stresses.push_back(tuple);
            //Increment group counter
            group_number++;
            //Reset average stress and datapoints number
            temp_stress_value = 0;
            datapoint_number = 0;
            }
        //Add the stress of an atom to a total
        temp_stress_value += System::surface_stresses[i][1];
        //Record the fact that some atom was found in the group
        datapoint_number++;
        prev_dist = System::surface_stresses[i][0];
        }
    //Record the lase group
    tuple.clear();
    tuple.push_back(group_number*binwidth);
    tuple.push_back(datapoint_number == 0 ? 0 : temp_stress_value/datapoint_number);
    averaged_stresses.push_back(tuple);

    //REcord the averaged stresses as a system member
    System::averaged_stresses = averaged_stresses;
    return averaged_stresses;
    }

//Calculate average stress of all atoms
float System::system_stress(){
    float total_stress = 0;
    for(Atom& atm : System::atoms){
        total_stress += atm.get_ave_stress();
        }
    return total_stress/System::atoms_number;
    }

//Isolate the surface to assess visually
void System::isolate_surface(array<string, 2> header, string filename){
    //Calculate new number of atoms
    string atoms_num = Helper::to_str((int)System::surface_atoms.size());
    //Recreate the appropriate header
    string contents = header[0] + atoms_num + "\n" + header[1];
    ofstream myfile;
    myfile.open(filename);
    myfile << contents;
    //Wrte down the surface atoms
    for(Atom& atm : System::surface_atoms){
        myfile << (atm.get_pure_line() + "\n");
        }
    myfile.close();
    }

//Find hydrogen species on the surface
map<string, float> System::get_surface_species(int htype){
    System::grid -> reset_grid(&(System::atoms), System::cutoffs);
    int H_ctr = 0;
    map<string, int> species_count;
    map<string, float> result;
    Molecule* ml;
    set<int> iter_exclude;
    Atom prev_atm = Atom();
    prev_atm.set_type(0);
    tuple<int, Molecule*> closestH;
    Molecule* truncated;
    for(Atom& atm : System::surface_atoms){
        iter_exclude.clear();
        ml = System::scan_molecule(prev_atm, atm, &iter_exclude);
        iter_exclude.clear();
        closestH = ml -> closest(htype);
        if(get<int>(closestH) < 3){
            truncated = get<Molecule*>(closestH) -> truncate(3);
            H_ctr++;
            species_count[truncated -> get_name()]++;
            for(int id : truncated -> get_ids()){
                System::exclude.insert(id);
                }
            }
        prev_atm = atm;
        }
    cout << "Number of Hydrogen species detected: " << H_ctr << endl;
    for(auto it = species_count.begin(); it != species_count.end(); ++it){
        result[it -> first] = (it -> second)/(float)System::surface_atoms.size();
        cout << it -> first << " " << result[it -> first] << endl;
        }
    return result;
    }


Molecule* System::scan_molecule(Atom prev_atm, Atom atm, set<int>* exclude){
//    cout << "Starting lookup for: " << atm.get_id() << "-----------------------" << endl;
    int atom_type = atm.get_type();
    Atom loop_atm;
    float dist;
    Molecule* res = new Molecule(atom_type, System::types[atom_type], atm.get_id());
    exclude -> insert(atm.get_id());
    vector<tuple<Atom, float>> neighbors = System::grid -> find_neighbors(&atm, *exclude, System::cutoffs);
    for(tuple<Atom, float>pair : neighbors){
        loop_atm = get<Atom>(pair);
        dist = get<float>(pair);
//        cout << "Current neighbor is: " << loop_atm.get_id() << endl;
        exclude -> insert(loop_atm.get_id());
        if(prev_atm.get_type() != 1 && loop_atm.get_type() != 3 && loop_atm.get_type() != 4){
            res -> add_bond(System::scan_molecule(atm, loop_atm, exclude), dist);
            }
        else{
            atom_type = loop_atm.get_type();
            if((atom_type != 3 && atom_type != 4) || loop_atm.Oid == atm.get_id()){
                res -> add_bond(atom_type, System::types[atom_type], loop_atm.get_id(), dist);
                }
            }
        }
//    cout << "Finished lookup for: " << atm.get_id() << "-----------------------" << endl;
    return res;
    }
