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
System::System(array<float, 3> box, ifstream& contents, int atoms_number, array<float, 3> center, array<float, 3> box_shift, int htype, int natype, float thickness){
    //Ionic settings
    System::htype = htype;
    System::natype = natype;
    System::surface_thickness = thickness;
    System::types[htype] = "H";
    System::types[natype] = "Na";
    System::radii_mapping[htype] = 0.53;
    System::radii_mapping[natype] = 2.27;
    System::cutoffs[2][htype] = 1.1;
    System::cutoffs[2][natype] = 2.7;
    System::cutoffs[htype][2] = 1.1;
    System::cutoffs[natype][2] = 2.7;

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
    int prev_num_surface = 0;
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
    int test_ctr = 0;
    //Loop while the grid volume is greater than minimum detected void
    while(!done){
        prev_num_surface = cur_num_surface;
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
            surface_atoms = System::grid -> get_surface(System::masks, System::surface_thickness);
            cur_num_surface = surface_atoms.size();
            //Increase the rastering definition
            num_splits *= 2;
            cout << "Current number of surface atoms: " << cur_num_surface << endl;
            cout << "Number of masks: " << masks -> get_member_num() << endl;
            cout << "-------------------------------\n";
            if(prev_num_surface && abs(cur_num_surface - prev_num_surface)/(float)prev_num_surface < 0.01){
//            cout << prev_num_surface << " " << cur_num_surface << " " << abs(cur_num_surface - prev_num_surface)/prev_num_surface << endl;
        done = true;
                }
            }
        else{done = true;}
        test_ctr++;
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
        atm -> set_potential(stof(split_line[11]));
        atm -> set_radius(System::radii_mapping[atm -> get_type()]);
        //Add the atom to the system
        System::atoms.push_back(*atm);
        //Record the modifier atoms in the modifiers vector
        if(atm -> get_type() != 1 && atm -> get_type() != 2 ){
            System::modifiers.push_back(*atm);
            }
        if(atm -> get_type() == 1){
            System::former_atoms.push_back(*atm);
            }
        if(atm -> get_type() == 2){
            System::oxygens.push_back(*atm);
            }
        if(atm -> get_type() == System::htype){
            System::hydrogens.push_back(*atm);
            }
        delete atm;
        }
    //A safe guard from incorrect atom number listed in .dump
    if(i < System::atoms_number){
        System::atoms_number = i;
        }
    }

//A function that calculates per atom potential energy as a function to the closest modifier
vector<vector<float>> System::calc_potentials(vector<Atom>& main_atoms, vector<Atom>& secondary_atoms){
    cout << "Beginning the potential energy calculations\n\n";
    tuple<Atom, float> closest;
    vector<vector<float>> stress_function;
    vector<float> tuple;
    AtomGrid* grid = new AtomGrid(System::box, 32, &secondary_atoms);
    //For each modifier atom
    for(Atom& atm : main_atoms){
        //Get closest atom
        closest = grid -> find_closest(&atm);
        tuple.clear();
        //Put the distance and stress into a resulting vector
        tuple.push_back(get<float>(closest));
        tuple.push_back(atm.get_potential());
        stress_function.push_back(tuple);
        }
    return stress_function;
    }

//A function that calculates per atom stress as a function to the closest modifier
vector<vector<float>> System::calc_stresses(vector<Atom>& main_atoms, vector<Atom>& secondary_atoms){
    cout << "Beginning the stress calculations\n\n";
    tuple<Atom, float> closest;
    vector<vector<float>> stress_function;
    vector<float> tuple;
    AtomGrid* grid = new AtomGrid(System::box, 32, &secondary_atoms);
    //For each modifier atom
    for(Atom& atm : main_atoms){
        //Get closest atom
        closest = grid -> find_closest(&atm);
        tuple.clear();
        //Put the distance and stress into a resulting vector
        tuple.push_back(get<float>(closest));
        tuple.push_back(atm.get_stress_comp(2));
        stress_function.push_back(tuple);
        }
    return stress_function;
    }

vector<vector<float>> System::average_property(vector<vector<float>>& stresses){
    /*
    sort(stresses.begin(), stresses.end(), [=](vector<float>& vect1, vector<float>& vect2){return vect1[0] < vect2[0];});
    float low_bound = (*stresses.begin())[0];
    float up_bound = (*stresses.end())[0];
    */
    float low_bound = (*min_element(stresses.begin(), stresses.end(), [=](vector<float>& vect1, vector<float>& vect2){return vect1[0] < vect2[0];}))[0];
    float up_bound = (*max_element(stresses.begin(), stresses.end(), [=](vector<float>& vect1, vector<float>& vect2){return vect1[0] < vect2[0];}))[0];
    cout << "Low: " << low_bound << " Up: " << up_bound << endl;
    return System::average_property(stresses, low_bound, up_bound);
    }

//Average the property by distance regions
vector<vector<float>> System::average_property(vector<vector<float>>& stresses, float low_bound, float up_bound){
    float binwidth = (up_bound - low_bound)/20;
    //Sort the stresses vector by distance to closest modifier
    sort(stresses.begin(), stresses.end(), [=](vector<float>& vect1, vector<float>& vect2){return vect1[0] < vect2[0];});
    float temp_stress_value = 0.0;
    int datapoint_number = 0, group_number = (int)(stresses[0][0]/binwidth);
    float prev_dist = -1.0;

    vector<vector<float>> averaged_stresses;
    vector<float> tuple;

    int stresses_vect_size = stresses.size();

    int i = 0;
	while(i < stresses_vect_size && stresses[i][0] < low_bound){
		i++;
		}

    //For each dist-stress pair
    for(; i < stresses_vect_size && stresses[i][0] < up_bound; i++){
        //If in the next distance region
        if(prev_dist < group_number*binwidth && stresses[i][0] > group_number*binwidth){
            tuple.clear();
            //Record average distance
            tuple.push_back(group_number*binwidth);
            //Record average stress
            tuple.push_back(datapoint_number == 0 ? 0 : temp_stress_value/datapoint_number);
            //Record the number of atoms contributing to this value
            tuple.push_back(datapoint_number);
            averaged_stresses.push_back(tuple);
            //Increment group counter
            group_number++;
            //Reset average stress and datapoints number
            temp_stress_value = 0;
            datapoint_number = 0;
            }
        //Add the stress of an atom to a total
        temp_stress_value += stresses[i][1];
        //Record the fact that some atom was found in the group
        datapoint_number++;
        prev_dist = stresses[i][0];
        }
    //Record the last group
    tuple.clear();
    tuple.push_back(group_number*binwidth);
    tuple.push_back(datapoint_number == 0 ? 0 : temp_stress_value/datapoint_number);
    tuple.push_back(datapoint_number);
    averaged_stresses.push_back(tuple);

    //REcord the averaged stresses as a system member
    System::averaged_stresses = averaged_stresses;
    return averaged_stresses;
    }

//Calculate average stress of all atoms
float System::system_stress(){
    float total_stress_x = 0;
    float total_stress_y = 0;
    float total_stress_z = 0;
    for(Atom& atm : System::atoms){
        total_stress_x += atm.get_stress_comp(0);
        total_stress_y += atm.get_stress_comp(1);
        total_stress_z += atm.get_stress_comp(2);
        }
    return pow(pow(total_stress_x, 2) + pow(total_stress_y, 2) + pow(total_stress_z, 2), 0.5)/(10*System::get_volume());
//    return total_stress_z/(10*System::get_volume());
    }

array<float, 2> System::average_potential(vector<Atom> atoms_studied){
    float potential = 0;
    float deviation_sum = 0;
    int numatoms = atoms_studied.size();
    for(Atom& atm : atoms_studied){
        potential += atm.get_potential();
        }
    float mean = potential/numatoms;
    for(Atom& atm : atoms_studied){
        deviation_sum += (float)pow(atm.get_potential() - mean, 2);
        }
    array<float, 2>result = {mean, (float)pow(deviation_sum/numatoms, 0.5)};
    return result;
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

/*
//Find hydrogen species on the surface
map<string, float> System::get_surface_species(){
    System::grid -> reset_grid(&(System::atoms), System::cutoffs);
    int mod_ctr = 0;
    map<string, int> species_count;
    map<string, float> result;
    Molecule* ml;
    set<int> iter_exclude;
    Atom prev_atm = Atom();
    prev_atm.set_type(0);
	int total_molecules = 0;
    tuple<int, Molecule*> closestH;
    tuple<int, Molecule*> closestNa;
    tuple<int, Molecule*> closestMod;
    Molecule* truncated;
    for(Atom& atm : System::surface_atoms){
        iter_exclude.clear();
        ml = System::scan_molecule(prev_atm, atm, &iter_exclude);
        iter_exclude.clear();
        closestH = ml -> closest(System::htype);
        closestNa = ml -> closest(System::natype);
		closestMod = get<int>(closestH) < get<int>(closestNa) ? closestH : closestNa;
        if(get<int>(closestMod) < 3){
            truncated = get<Molecule*>(closestMod) -> truncate(3);
            mod_ctr++;
            species_count[truncated -> get_name()]++;
            for(int id : truncated -> get_ids()){
                System::exclude.insert(id);
                if(truncated -> get_name() == "O(1)H(1)"){
                    cout << id << endl;
                    }
                }
            }
        prev_atm = atm;
		total_molecules++;
        }
    cout << "Number of Modifier species detected: " << mod_ctr << endl;
    for(auto it = species_count.begin(); it != species_count.end(); ++it){
        result[it -> first] = (it -> second)/(float)total_molecules;
        cout << it -> first << " " << result[it -> first] << endl;
        }
    return result;
    }
*/

//Oxygen neighbor analysis
map<string, float> System::get_surface_species(){
    System::grid -> reset_grid(&(System::atoms), System::cutoffs);
    map<string, int> species_count;
    map<string, float> result;
    set<int> exclude;
    map<int, int> oxygen_neighbors;
    string name;
    vector<Atom> surface_oxygens = System::get_species(System::surface_atoms, 2);
    vector<tuple<Atom, float>> neighbors;
    int ox_id, si_id;
    int neigh_type;
    tuple<Atom, float> closest;
    for(Atom& atm : surface_oxygens){
        exclude.clear();
        name = "";
        ox_id = atm.get_id();
        oxygen_neighbors.clear();
        neighbors = System::grid -> find_neighbors(&atm, exclude, System::cutoffs);
        for(tuple<Atom, float>pair : neighbors){
            neigh_type = get<Atom>(pair).get_type();
            if(neigh_type == 1){si_id = get<Atom>(pair).get_id();}
            if((neigh_type != System::htype && neigh_type != System::natype) || get<Atom>(pair).Oid == ox_id){
                oxygen_neighbors[neigh_type]++;
                }
            }
        for(auto it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); ++it){
            name += (System::types[it -> first] + "(" + Helper::to_str(it -> second) + ")");
            }
        if(name == ""){
            name = "(-)";
            }
        if(name == "Si(1)"){
            exclude.insert(si_id);
            closest = System::grid -> find_closest(&atm, exclude);
            if(get<float>(closest) <= 1.3*(System::radii_mapping[2] + System::radii_mapping[get<Atom>(closest).get_type()])){
                species_count["NBO - " + System::types[get<Atom>(closest).get_type()]]++;
                }
            }
        species_count[name]++;
        }
    int total_oxygens = surface_oxygens.size();
    for(auto it = species_count.begin(); it != species_count.end(); ++it){
        result[it -> first] = (it -> second)/(float)total_oxygens;
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
        atom_type = loop_atm.get_type();
        //Continue if did not hit modifier (+1 charge) and if silicon is not the previous atom.
        if(prev_atm.get_type() != 1 && atom_type != System::htype && atom_type != System::natype){
            res -> add_bond(System::scan_molecule(atm, loop_atm, exclude), dist);
            }
        //Add terminating bond to silicon neighbor or correct oxygen neighbor
        else{
            if(atom_type != 3 && atom_type != 4){
                res -> add_bond(atom_type, System::types[atom_type], loop_atm.get_id(), dist);
                }
            }
        }
//    cout << "Finished lookup for: " << atm.get_id() << "-----------------------" << endl;
    return res;
    }

vector<Atom> System::get_hydrogens(){
    return System::hydrogens;
    }

vector<Atom> System::get_modifiers(){
    return System::modifiers;
    }

vector<Atom> System::get_formers(){
    return System::former_atoms;
    }

vector<Atom> System::get_oxygens(){
    return System::oxygens;
    }

float System::get_volume(){
    return System::box[0]*System::box[1]*System::box[2];
    }

int System::count_species(vector<Atom> sample, int type){
    int ctr = 0;
    for(Atom& atm : sample){
        if(atm.get_type() == type){
            ctr++;
            }
        }
    return ctr;
    }

float System::get_total_comp(int index){
    float total = 0;
    for(Atom& atm : System::atoms){
        total += atm.get_stress_comp(index);
        }
    return total/(10*System::get_volume());
    }

vector<vector<float>> System::get_peratom_stress(int index){
    vector<vector<float>> result;
    vector<float> atom_stress;
    for(Atom& atm : System::atoms){
        atom_stress.clear();
        atom_stress.push_back((float)(atm.get_id()));
        atom_stress.push_back((float)(atm.get_stress_comp(index)));
        result.push_back(atom_stress);
        }
    return result;
    }

vector<Atom> System::get_species(vector<Atom> studied_atoms, int type){
	vector<Atom> result;
	for(Atom& atm : studied_atoms){
		if(atm.get_type() == type){
			result.push_back(atm);
			}
		}
	return result;
	}

float System::get_total_potential(vector<Atom> studied_atoms){
	float potential = 0;
	for(Atom& atm : studied_atoms){
		potential += atm.get_potential();
		}
	return potential;
	}

map<int, float> System::get_Qunits(){
    cout << "Q-unit calculations\n";
    System::grid -> reset_grid(&(System::atoms), System::cutoffs);
    map<int, int> qunits;
    vector<tuple<Atom, float>> si_neighbors;
    vector<tuple<Atom, float>> o_neighbors;
    Atom si_neigh;
    int q_count;
    int num_si = System::former_atoms.size();
    set<int> exclude;
    int test_ctr = 0;
    float percent_done;
    printf("\e[?25l");
    for(Atom &atm : System::former_atoms){
//        cout << "\r" << "Percent q-unit scan complete: ";
//        printf("%.2f", 100.0*test_ctr/num_si);
        exclude.clear();
        exclude.insert(atm.get_id());
        q_count = 0;
        si_neighbors = System::grid -> find_neighbors(&atm, exclude, System::cutoffs);
        for(tuple<Atom, float> &pair : si_neighbors){
            si_neigh = get<Atom>(pair);
            if(si_neigh.get_type() == 2){
                o_neighbors = System::grid -> find_neighbors(&si_neigh, exclude, System::cutoffs);
                for(tuple<Atom, float> &o_pair : o_neighbors){
                    if(get<Atom>(o_pair).get_type() == 1){
                        q_count++;
                        }
                    }
                }
            }
        test_ctr++;
        qunits[q_count]++;
        }
    cout << endl;
    printf("\e[?25h");
    map<int, float> q_conc;
    for(auto it = qunits.begin(); it != qunits.end(); ++it){
        q_conc[it -> first] = (float)(it -> second)/num_si;
        }
    return q_conc;
    }
