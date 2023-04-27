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

string get_box(ifstream& contents, array<float, 3>& box, array<float, 3>& box_shift, array<float, 3>& center){
    string combined_lines = "";
    string line;
    vector<string> split_line;
    for(int i = 0; i < 3; i++){
        getline(contents, line);
        combined_lines += (line + "\n");
        split_line = Helper::split(line, " ");
        box_shift[i] = stof(split_line[0]);
        box[i] = stof(split_line[1]) - box_shift[i];
        center[i] = box[i]/2;
        }
    return combined_lines;
    }
vector<float> analysis(string destination, string filename, bool iso_surface){
    ifstream contents(destination + filename);
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
    while(getline(contents, line)){
        header[header_ctr] += (line + "\n");
        if(Helper::string_contains(line, "ITEM: NUMBER OF ATOMS")){
            getline(contents, line);
            atoms_number = stoi(line);
            header_ctr++;
            }
        if(Helper::string_contains(line, "ITEM: BOX BOUNDS")){
            header[header_ctr] += get_box(contents, box, box_shift, center);
            }
        if(Helper::string_contains(line, "ITEM: ATOMS")){
            System* atom_system = new System(box, contents, atoms_number, center, box_shift);
            total_stresses = atom_system -> calc_stresses();
            Helper::vector2d_csv(destination + "total/" + filename, "Distance to closest modifier, Stress", total_stresses);
            average_stresses = atom_system -> average_stresses(binwidth);
            Helper::vector2d_csv(destination + "averaged/" + filename, "Bin span, Stress", average_stresses);
            if(iso_surface){
                atom_system -> isolate_surface(header, "surface." + filename);
                }
            result.push_back(box[2]);
            result.push_back(atom_system -> system_stress());
            }
        }
    return result;
    }

int main(int argc, char** argv){
    string input(argv[1]);
    fs::path cwd = fs::current_path();
    if(Helper::string_contains(input, "*")){
        vector<vector<float>> stresses;
        if(Helper::string_contains(input, "/")){
            vector<string> split_name = Helper::split(input, "/");
            string pattern = split_name.back();
            string file_location = "";
            for(string& dir_name : split_name){
                file_location += (dir_name + "/");
                }
            string destination = "stress_analysis";
            fs::create_directory(cwd/destination/"total");
            fs::create_directory(cwd/destination/"averaged");
            for(string& filename : Helper::files_by_pattern(file_location, pattern)){
                stresses.push_back(analysis(destination + "/", filename, false));
                }
            }
        }else{
            fs::create_directory(cwd/"total");
            fs::create_directory(cwd/"averaged");
            analysis("./", input, true);
            }
    return 0;
    }
