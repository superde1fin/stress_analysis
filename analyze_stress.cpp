#include <iostream>
#include <string>
#include <filesystem>
#include <array>
#include <vector>

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
vector<float> analysis(string file_location, string destination, string filename, bool iso_surface, float void_volume){
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
    array<string, 2> header = {"", ""};
    int header_ctr = 0;
    vector<Atom> surface;
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
            System* atom_system = new System(box, contents, atoms_number, center, box_shift);
            //Run the surface detection
            surface = atom_system -> detect_surface(void_volume);
            //If specified, create the surface file
            if(iso_surface){
                atom_system -> isolate_surface(header, destination + "surfaces/surface." + filename);
                }
            //Run the stress calculator
            total_stresses = atom_system -> calc_stresses();
            //Record the atom stresses into a csv
            Helper::vector2d_csv(destination + "total/" + filename, "Distance to closest modifier, Stress", total_stresses);
            //Condense the data and average stress based on distance ranges
            average_stresses = atom_system -> average_stresses(binwidth);
            Helper::vector2d_csv(destination + "averaged/" + filename, "Bin span, Stress", average_stresses);
            //Record the box dimension and average stress in the result to create the stress-strain curve
            result.push_back(box[2]);
            result.push_back(atom_system -> system_stress());
            }
        }
    return result;
    }

int main(int argc, char** argv){
    //Read file name
    string input(argv[1]);
    //Read minimum void space
    string void_volume(argv[2]);
    //Create the file system object
    fs::path cwd = fs::current_path();
    //Check if the inputted name is a pattern
    if(Helper::string_contains(input, "*")){
        vector<vector<float>> stresses;
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
            fs::create_directory(cwd/destination/"total");
            fs::create_directory(cwd/destination/"averaged");
            fs::create_directory(cwd/destination/"surfaces");
            //For each file that matches the pattern perform the stress calculations
            for(string& filename : Helper::files_by_pattern(file_location, pattern, true)){
                cout << "Starting the stress analysis for: " << filename << endl << endl;
                //Store the stress and strain in a csv file
                stresses.push_back(analysis(file_location, destination + "/", filename, true, stof(void_volume)));
                Helper::vector2d_csv(destination + "/stress_strain", "Box z, Ave Stress", stresses);
                }
            }
        }else{
            //One file case
            fs::create_directory(cwd/"total");
            fs::create_directory(cwd/"averaged");
            fs::create_directory(cwd/"surfaces");
            cout << "Starting the stress analysis for: " << input << endl << endl;
            analysis("./", "./", input, true, stof(void_volume));
            }
    return 0;
    }
