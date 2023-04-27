#include <array>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <limits>

#include "System.hpp"
#include "Atom.hpp"
#include "Ellipse.hpp"
#include "Sector.hpp"
#include "Helper.hpp"

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

vector<Atom> System::detect_surface(){
    vector<Atom> surface_atoms;
    float step = System::min_radius;
    bool success, done;
    int i, atom_list_len;
    float dist_init, elliptical_dist;
    vector<vector<Atom>> atom_slabs = System::split_atoms(step);
    for(vector<Atom>& atom_list : atom_slabs){
        sort(atom_list.begin(), atom_list.end(), [=](Atom& atm1,Atom& atm2){
                return pow(pow(atm1.get_y() - System::center[1], 2) + 25*pow(atm1.get_z() - System::center[2], 2), 0.5) < pow(pow(atm2.get_y() - System::center[1], 2) + 25*pow(atm2.get_z() - System::center[2], 2), 0.5);
                });
        atom_list_len = atom_list.size();
        Ellipse* crcl;
        Sector* seg;
        i = 0;
        dist_init = pow(pow(atom_list[i].get_y() - System::center[1], 2) + 25*pow(atom_list[i].get_z() - System::center[2], 2), 0.5)/5;
        elliptical_dist = dist_init;
        done = false;
        Sector::reset();
        while(!done){
//            cout << dist_init << " " << elliptical_dist << " Center: " << System::center[1] << " " << System::center[2] << endl << atom_list[i].get_y() << " " << atom_list[i].get_z() << endl;
            if(atom_list_len && i < atom_list_len && elliptical_dist < dist_init + 5*step){
//            if(atom_list_len && i < 3){
//            if(atom_list_len && i < atom_list_len){
                crcl = new Ellipse(System::center, step, atom_list[i]);
                seg = new Sector(*crcl, atom_list[i], &success);
                if(success){
                    surface_atoms.push_back(atom_list[i]);
                    elliptical_dist = crcl -> get_minor();
                    if(Sector::get_area(*crcl)/(crcl -> get_area()) >= 0.9){done = true;}
                    }
                }else{done = true;}

            i++;
            }
//            break;
        }
    System::surface_atoms = surface_atoms;
    return surface_atoms;
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
        atm -> set_position(array<float, 3>{stof(split_line[2]) - System::box_shift[0], stof(split_line[3]) - System::box_shift[1], stof(split_line[4])});
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


vector<vector<Atom>> System::split_atoms(float step){
    int num_regions = (int)(System::box[0]/step) + 1;
    vector<vector<Atom>> slabs(num_regions);
    for(int i = 0; i < System::atoms_number; i++){
        slabs[(int)(fmod(System::atoms[i].get_x(), System::box[0])/step)].push_back(System::atoms[i]);
        }
    return slabs;
    }

vector<vector<float>> System::calc_stresses(){
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
    for(Atom& atm : System::surface_atoms){
        contents += (atm.get_pure_line() + "\n");
        }
    ofstream myfile;
    myfile.open(filename);
    myfile << contents;
    myfile.close();
    }
