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
tuple<vector<float>, map<string, float>> analysis(string file_location, string destination, string filename, bool iso_surface, float void_volume, int htype, int natype){
    ifstream contents(file_location + filename);
    string line;
    int atoms_number;
    float binwidth = 0.2;
    array<float, 3> box;
    array<float, 3> box_shift;
    array<float, 3> center;
    vector<vector<float>> total_stresses;
    vector<vector<float>> average_stresses;
    vector<vector<float>> peratom_stresses;
    vector<float> result;
    map<string, float> species;
    array<string, 2> header = {"", ""};
    int header_ctr = 0;
    vector<Atom> modifiers;
    vector<Atom> hydrogens;
    vector<Atom> surface; 
    vector<Atom> bulk; 
    array<float, 2> potential_result;
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
            System* atom_system = new System(box, contents, atoms_number, center, box_shift, htype, natype);
            //Run the surface detection
            surface = atom_system -> detect_surface(void_volume);
            //If specified, create the surface file
            if(iso_surface){
                atom_system -> isolate_surface(header, destination + "surfaces/surface." + filename);
                }
//-----------------SURFACE
            //MODIFIERS BASED STRESS
            modifiers = atom_system -> get_modifiers();
            //Run the stress calculator based on modifiers
            total_stresses = atom_system -> calc_stresses(surface, modifiers);
            //Record the atom stresses into a csv
            Helper::vector2d_csv(destination + "surface/modifiers/total/" + filename, "Distance to closest modifier, Potential Energy", total_stresses);
            //Condense the data and average stress based on distance ranges
            average_stresses = atom_system -> average_stresses(total_stresses, 2.0, 4.0);
            Helper::vector2d_csv(destination + "surface/modifiers/averaged/" + filename, "Bin span, Potential Energy, Atom count", average_stresses);
            result.push_back(box[2]);
            potential_result = atom_system -> average_potential(atom_system -> get_formers());
            result.push_back(potential_result[0]);
            result.push_back(potential_result[1]);
            potential_result = atom_system -> average_potential(atom_system -> get_oxygens());
            result.push_back(potential_result[0]);
            result.push_back(potential_result[1]);
            result.push_back((float)(atom_system -> count_species(surface, htype)));
            result.push_back((float)(atom_system -> count_species(surface, natype)));
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
            fs::create_directory(cwd/destination/"bulk");
            fs::create_directory(cwd/destination/"surface");
//            fs::create_directory(cwd/destination/"XX_stresses");
//            fs::create_directory(cwd/destination/"YY_stresses");
//            fs::create_directory(cwd/destination/"ZZ_stresses");
//            fs::create_directory(cwd/destination/"bulk/modifiers");
//            fs::create_directory(cwd/destination/"bulk/hydrogens");
//            fs::create_directory(cwd/destination/"bulk/modifiers/total");
//            fs::create_directory(cwd/destination/"bulk/modifiers/averaged");
//            fs::create_directory(cwd/destination/"bulk/hydrogens/total");
//            fs::create_directory(cwd/destination/"bulk/hydrogens/averaged");
            fs::create_directory(cwd/destination/"surface/modifiers");
//            fs::create_directory(cwd/destination/"surface/hydrogens");
            fs::create_directory(cwd/destination/"surface/modifiers/total");
            fs::create_directory(cwd/destination/"surface/modifiers/averaged");
//            fs::create_directory(cwd/destination/"surface/hydrogens/total");
//            fs::create_directory(cwd/destination/"surface/hydrogens/averaged");
//            fs::create_directory(cwd/destination/"surfaces");
            //For each file that matches the pattern perform the stress calculations
            vector<string> pattern_fillers;
            for(string& filename : Helper::files_by_pattern(file_location, pattern, &pattern_fillers, true)){
                cout << "Starting the stress analysis for: " << filename << endl << endl;
                analysis_result = analysis(file_location, destination + "/", filename, false, stof(void_volume), stoi(htype), stoi(natype));
                //Store the stress and strain in a csv file
                system_output.push_back(get<vector<float>>(analysis_result));
                surface_species.push_back(get<map<string, float>>(analysis_result));
                Helper::vector2d_csv(destination + "/system_output", "Timestep, Box z, Average Potential Energy on Surface Si, Standard Deviation, Average Potential Energy on Surface O, Standard Deviation, H num, Na num", system_output, pattern_fillers);
//                Helper::vector_of_maps2csv(destination + "/species", surface_species, pattern_fillers);
                }
            }
        }else{
            //One file case
            fs::create_directory(cwd/"bulk");
            fs::create_directory(cwd/"surface");
            fs::create_directory(cwd/"XX_stresses");
            fs::create_directory(cwd/"YY_stresses");
            fs::create_directory(cwd/"ZZ_stresses");
            fs::create_directory(cwd/"bulk/modifiers");
            fs::create_directory(cwd/"bulk/hydrogens");
            fs::create_directory(cwd/"bulk/modifiers/total");
            fs::create_directory(cwd/"bulk/modifiers/averaged");
            fs::create_directory(cwd/"bulk/hydrogens/total");
            fs::create_directory(cwd/"bulk/hydrogens/averaged");
            fs::create_directory(cwd/"surface/modifiers");
            fs::create_directory(cwd/"surface/hydrogens");
            fs::create_directory(cwd/"surface/modifiers/total");
            fs::create_directory(cwd/"surface/modifiers/averaged");
            fs::create_directory(cwd/"surface/hydrogens/total");
            fs::create_directory(cwd/"surface/hydrogens/averaged");
            fs::create_directory(cwd/"surfaces");
            cout << "Starting the stress analysis for: " << input << endl << endl;
            analysis("./", "./", input, true, stof(void_volume), stoi(htype), stoi(natype));
            }
    return 0;
    }
