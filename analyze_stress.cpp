#include <iostream>
#include <string>
#include <filesystem>
#include <array>
#include <vector>
#include <tuple>

#include "classes/System.hpp"
#include "classes/Atom.hpp"
#include "classes/Helper.hpp"

using namespace std;
namespace fs = std::filesystem;

//This function parses the part of a .dump file and returns the box dimensions
string get_box(ifstream& contents, array<float, 3>& box, array<float, 3>& box_shift, array<float, 3>& center){
    string combined_lines = "";
    string line;
    vector<string> split_line;
    for(int i = 0; i < 3; i++){
        //read line
        getline(contents, line);
        //record for future surface recreation
        combined_lines += (line + "\n");
        //split line by spaces
        split_line = Helper::split(line, " ");
        //record box shift (lammps generated shift of coordinates)
        box_shift[i] = stof(split_line[0]);
        //record box dimensions from 0
        box[i] = stof(split_line[1]) - box_shift[i];
        //record center coordinate
        center[i] = box[i]/2;
        }
    return combined_lines;
    }
tuple<vector<float>, map<string, float>> analysis(string file_location, string destination, string filename, bool iso_surface, float void_volume, int htype, int natype, float thickness){
    ifstream contents(file_location + filename);
    string line;
    int atoms_number;
    float binwidth = 0.2;
    array<float, 3> box;
    array<float, 3> box_shift;
    array<float, 3> center;
    vector<vector<float>> total_stresses;
    vector<vector<float>> average_stresses;
    vector<float> result;
    map<string, float> species;
    array<string, 2> header = {"", ""};
    int header_ctr = 0;
    vector<Atom> modifiers;
    vector<Atom> atomic_species;
    vector<Atom> surface; 
    array<float, 2> potential_result;
    map<int, float> qunits;
    //Scan through lines
    while(getline(contents, line)){
        //Record header for surface recreation
        header[header_ctr] += (line + "\n");
        if(Helper::string_contains(line, "ITEM: NUMBER OF ATOMS")){
            getline(contents, line);
            //Record the number of atoms in the system
            atoms_number = stoi(line);
            header_ctr++;
            }
        if(Helper::string_contains(line, "ITEM: BOX BOUNDS")){
            //Call the box dimensions parser
            header[header_ctr] += get_box(contents, box, box_shift, center);
            }
        if(Helper::string_contains(line, "ITEM: ATOMS")){
            //Create a system object
            System* atom_system = new System(box, contents, atoms_number, center, box_shift, htype, natype, thickness);
            //Run the surface detection
            surface = atom_system -> detect_surface(void_volume);
            //If specified, create the surface file
            if(iso_surface){
                atom_system -> isolate_surface(header, destination + "surfaces/surface." + filename);
                }

            result.push_back(box[2]);
//----------Stresses on surface silicons
            modifiers = atom_system -> get_modifiers();
			atomic_species = atom_system -> get_species(surface, 1);
            //Run the stress calculator based on modifiers
            total_stresses = atom_system -> calc_stresses(atomic_species, modifiers);
            //Record the atom stresses into a csv
            Helper::vector2d_csv(destination + "stresses/silicon/total/" + filename, "Distance to closest modifier, ZZ stress", total_stresses);
            //Condense the data and average stress based on distance ranges
            average_stresses = atom_system -> average_property(total_stresses, 2.0, 4.0);
            Helper::vector2d_csv(destination + "stresses/silicon/averaged/" + filename, "Bin span, ZZ stress, Atom count", average_stresses);
//----------Potential on surface silicons
            //Run the stress calculator based on modifiers
            total_stresses = atom_system -> calc_potentials(atomic_species, modifiers);
            //Record the atom stresses into a csv
            Helper::vector2d_csv(destination + "potentials/silicon/total/" + filename, "Distance to closest modifier, Potential Energy", total_stresses);
            //Condense the data and average stress based on distance ranges
            average_stresses = atom_system -> average_property(total_stresses, 2.0, 4.0);
            Helper::vector2d_csv(destination + "potentials/silicon/averaged/" + filename, "Bin span, PotEng, Atom count", average_stresses);

//----------Stresses on surface oxygens
			atomic_species = atom_system -> get_species(surface, 2);
            //Run the stress calculator based on modifiers
            total_stresses = atom_system -> calc_stresses(atomic_species, modifiers);
            //Record the atom stresses into a csv
            Helper::vector2d_csv(destination + "stresses/oxygen/total/" + filename, "Distance to closest modifier, ZZ stress", total_stresses);
            //Condense the data and average stress based on distance ranges
            average_stresses = atom_system -> average_property(total_stresses, 1.0, 7.0);
            Helper::vector2d_csv(destination + "stresses/oxygen/averaged/" + filename, "Bin span, ZZ stress, Atom count", average_stresses);
//----------Potential on surface oxygens
            //Run the stress calculator based on modifiers
            total_stresses = atom_system -> calc_potentials(atomic_species, modifiers);
            //Record the atom stresses into a csv
            Helper::vector2d_csv(destination + "potentials/oxygen/total/" + filename, "Distance to closest modifier, Potential Energy", total_stresses);
            //Condense the data and average stress based on distance ranges
            average_stresses = atom_system -> average_property(total_stresses, 1.0, 7.0);
            Helper::vector2d_csv(destination + "potentials/oxygen/averaged/" + filename, "Bin span, PotEng, Atom count", average_stresses);

            result.push_back(atom_system -> get_total_comp(2));
            result.push_back((float)atom_system -> count_species(surface, 1));
            result.push_back((float)atom_system -> count_species(surface, 2));
            result.push_back((float)atom_system -> count_species(surface, htype));
            result.push_back((float)atom_system -> count_species(surface, natype));

            qunits = atom_system -> get_Qunits();
			if(qunits.count(1)){
            	result.push_back(qunits[1]);
				}
			else{
            	result.push_back(0);
				}
			if(qunits.count(2)){
            	result.push_back(qunits[2]);
				}
			else{
            	result.push_back(0);
				}
			if(qunits.count(3)){
            	result.push_back(qunits[3]);
				}
			else{
            	result.push_back(0);
				}
			if(qunits.count(4)){
            	result.push_back(qunits[4]);
				}
			else{
            	result.push_back(0);
				}

			//Get a map of molecular species on the surface
			species = atom_system -> get_surface_species();
            }
        }
    return make_tuple(result, species);
    }

int main(int argc, char** argv){
    //Read file name
    string input(argv[1]);
    //Read minimum void space
    string void_volume(argv[2]);
    //Record hydrogen type
    string htype(argv[3]);
    //Record sodium type
    string natype(argv[4]);
    //Get surface thickness
    string thickness;
    if(argc > 5){
        thickness = string(argv[5]);
        }
    else{
        thickness = "0";
        }
    //Create the file system object
    fs::path cwd = fs::current_path();
    tuple<vector<float>, map<string, float>> analysis_result;
    vector<map<string, float>> surface_species;
    vector<vector<float>> system_output;
    //Check if the inputted name is a pattern
    if(Helper::string_contains(input, "*")){
        //Check whether the files are in a subdirectory
        if(Helper::string_contains(input, "/")){
            //Split name by slashes
            vector<string> split_name = Helper::split(input, "/");
            //Get the actual file pattern
            string pattern = split_name.back();
            split_name.pop_back();
            string file_location = "";
            //Reconstruct the directory
            for(string& dir_name : split_name){
                file_location += (dir_name + "/");
                }
            //Create the subdirectories to store the results
            string destination = "analysis";
            fs::create_directory(cwd/destination);
            fs::create_directory(cwd/destination/"potentials");
            fs::create_directory(cwd/destination/"stresses");
            fs::create_directory(cwd/destination/"surfaces");
            fs::create_directory(cwd/destination/"potentials/silicon");
            fs::create_directory(cwd/destination/"stresses/silicon");
            fs::create_directory(cwd/destination/"potentials/oxygen");
            fs::create_directory(cwd/destination/"stresses/oxygen");
            fs::create_directory(cwd/destination/"potentials/silicon/averaged");
            fs::create_directory(cwd/destination/"potentials/silicon/total");
            fs::create_directory(cwd/destination/"stresses/silicon/averaged");
            fs::create_directory(cwd/destination/"stresses/silicon/total");
            fs::create_directory(cwd/destination/"potentials/oxygen/averaged");
            fs::create_directory(cwd/destination/"potentials/oxygen/total");
            fs::create_directory(cwd/destination/"stresses/oxygen/averaged");
            fs::create_directory(cwd/destination/"stresses/oxygen/total");
            //For each file that matches the pattern perform the stress calculations
            vector<string> pattern_fillers;
            for(string& filename : Helper::files_by_pattern(file_location, pattern, &pattern_fillers, true)){
                cout << "Starting the stress analysis for: " << filename << endl << endl;
                analysis_result = analysis(file_location, destination + "/", filename, true, stof(void_volume), stoi(htype), stoi(natype), stof(thickness));
                //Store the stress and strain in a csv file
                system_output.push_back(get<vector<float>>(analysis_result));
                surface_species.push_back(get<map<string, float>>(analysis_result));
                Helper::vector2d_csv(destination + "/system_output", "Timestep, Box z, Stress ZZ, Si num, O num, H num, Na num, Q1, Q2, Q3, Q4", system_output, pattern_fillers);
                Helper::vector_of_maps2csv(destination + "/species", surface_species, pattern_fillers);
                }
            }
        }else{
            //One file case
            string destination = "analysis";
            fs::create_directory(cwd/destination);
            fs::create_directory(cwd/destination/"potentials");
            fs::create_directory(cwd/destination/"stresses");
            fs::create_directory(cwd/destination/"surfaces");
            fs::create_directory(cwd/destination/"potentials/silicon");
            fs::create_directory(cwd/destination/"stresses/silicon");
            fs::create_directory(cwd/destination/"potentials/oxygen");
            fs::create_directory(cwd/destination/"stresses/oxygen");
            fs::create_directory(cwd/destination/"potentials/silicon/averaged");
            fs::create_directory(cwd/destination/"potentials/silicon/total");
            fs::create_directory(cwd/destination/"stresses/silicon/averaged");
            fs::create_directory(cwd/destination/"stresses/silicon/total");
            fs::create_directory(cwd/destination/"potentials/oxygen/averaged");
            fs::create_directory(cwd/destination/"potentials/oxygen/total");
            fs::create_directory(cwd/destination/"stresses/oxygen/averaged");
            fs::create_directory(cwd/destination/"stresses/oxygen/total");
            cout << "Starting the stress analysis for: " << input << endl << endl;
            analysis("./", destination + "/", input, true, stof(void_volume), stoi(htype), stoi(natype), stof(thickness));
            }
    //analysis(file_location, destination + "/", filename, true, stof(void_volume), stoi(htype), stoi(natype), stof(thickness));
    return 0;
    }
