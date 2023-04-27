#include <iostream>
#include <string>
#include <filesystem>

#include "classes/System.hpp"
#include "classes/Atom.hpp"
#include "classes/Helper.hpp"

using namespace std;
namespace fs = std::filesystem;

vector<float> analysis(string destination, string filename, bool iso_surface){
    ifstream contents(destination + filename);
    string line;
    int atoms_number;
    while(getline(contents, line)){
        if(Helper::string_contains(line, "ITEM: NUMBER OF ATOMS")){
            getline(contents, line);
            atoms_number = stoi(line);
            cout << atoms_number << endl;
            }
        }
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
            fs::create_directory(cwd/"stress_analysis/total");
            fs::create_directory(cwd/"stress_analysis/averaged");
            for(string& filename : Helper::files_by_pattern(file_location, pattern)){
                stresses.push_back(analysis("stress_analysis", filename, true));
                }
            }
        }else{
            fs::create_directory(cwd/"total");
            fs::create_directory(cwd/"averaged");
            analysis("./", input, false);
            }
    }
