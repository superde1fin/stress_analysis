#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <iostream>
#include <tuple>
#include <map>

#include "Helper.hpp"
#include "Atom.hpp"

using namespace std;
namespace fs = filesystem;

int Helper::key2int(array<int, 3> key){
    return pow(2, key[0])*pow(3, key[1])*pow(5, key[2]);
    }

void Helper::remove_dupl(vector<Atom>& atoms){
    map<int, int> atom_ctrs;
    vector<Atom> result;
    for(auto it = atoms.begin(); it != atoms.end(); it++){
        if(!atom_ctrs[(*it).get_id()]){
            result.push_back(*it);
            }
        atom_ctrs[(*it).get_id()]++;
        }
    atoms = result;
    }

vector<string> Helper::split(string str, string delim){
    vector<string> result;
    int str_size = str.length();
    int delim_index = 0;
    int delim_size = delim.length();
    string part = "";
    for(int i = 0; i < str_size; i++){
        if(str[i] == delim[delim_index]){
            if(delim_index == delim_size - 1){
                result.push_back(part);
                part = "";
                delim_index = 0;
                }else{delim_index++;}
            }else{
                part += str[i];
                delim_index = 0;
                }
        }
    result.push_back(part);
    return result;
    }

bool Helper::string_contains(string to_check, string substr){
    int string_length = to_check.length();
    string current_substr("");
    int substr_ctr = 0;
    for(int i = 0; i < string_length; i++){
        if(to_check[i] == substr[substr_ctr]){
            substr_ctr++;
            current_substr += to_check[i];
            if(current_substr == substr){return true;}
            }else{
                substr_ctr = 0;
                current_substr = "";
                }
        }
    return false;
    }

float Helper::dist(array<float, 3>atom1, array<float, 3> atom2, array<float, 3> box){
    float x, y, z; x = abs(atom1[0] - atom2[0]);
    x = (x > box[0]/2) ? (box[0] - x) : x;
    
    y = abs(atom1[1] - atom2[1]);
    y = (y > box[1]/2) ? (box[1] - y) : y;
    
    z = abs(atom1[2] - atom2[2]);
    z = (z > box[2]/2) ? (box[2] - z) : z;
    
    return pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5);
    }

tuple<Atom, float> Helper::find_closest(Atom atm, vector<Atom> to_compare, array<float, 3> box){
    float min = numeric_limits<float>::infinity();
    float cur_dist;
    Atom key;
    for(Atom& atmi : to_compare){
        cur_dist = dist(atm.get_position(), atmi.get_position(), box);
        if(cur_dist < min && cur_dist != 0){
            min = cur_dist; 
            key = atmi;
            }
        }
    return make_tuple(key, min);
    }

void Helper::vector2d_csv(string name, string first_line, vector<vector<float>> vect){
    char comma;
    float word;
    int line_length;
    vector<float> line;

    first_line += '\n';
    for (int j = 0; j < (int)vect.size(); j++){
        line = vect[j];
        comma = ',';
        line_length = (int)line.size();
        for (int i = 0; i < line_length; i++){
            word = line[i];
            if (i == line_length - 1){comma = '\n';}
            first_line += to_str(word) + comma;
            }
        }
    ofstream myfile;
    myfile.open(name + ".csv");
    myfile << first_line;
    myfile.close();
    }


bool Helper::fits_pattern(string to_check, string pattern){
    int string_length = to_check.length();
    int pattern_length = pattern.length();
    int i;
    for(i = 0; i < string_length && pattern[i] != '*'; i++){
        if(to_check[i] != pattern[i]){return false;}
        }
    if(pattern[i] != '*'){return false;}
    for(i = 1; i <= pattern_length && pattern[pattern_length - i] != '*'; i++){
        if(to_check[string_length - i] != pattern[pattern_length - i]){return false;}
        }
    if(pattern[pattern_length - i] != '*'){return false;}
    return true;
    }

bool Helper::fits_pattern(string to_check, string pattern, string* filler){
    int string_length = to_check.length();
    int pattern_length = pattern.length();
    int i, j;
    for(i = 0; i < string_length && pattern[i] != '*'; i++){
        if(to_check[i] != pattern[i]){return false;}
        }
    if(pattern[i] != '*'){return false;}
    for(j = 1; j <= pattern_length && pattern[pattern_length - j] != '*'; j++){
        if(to_check[string_length - j] != pattern[pattern_length - j]){return false;}
        }
    if(pattern[pattern_length - j] != '*'){return false;}
    string within("");
    for(int k = i; k <= string_length - j; k++){
        within += to_check[k];
        }
    *filler = within;
    return true;
    }

vector<string> Helper::files_by_pattern(string lookup_dir, string pattern){
    vector<string> found_files;
    string filename;
    string within;
    for (const auto & entry : fs::directory_iterator(lookup_dir)){
        filename = split(entry.path().string(), "/").back();
        if(fits_pattern(filename, pattern, &within)){
            found_files.push_back(filename);
            }
        }
    return found_files;
    }

vector<string> Helper::files_by_pattern(string lookup_dir, string pattern, bool sort_results = false){
    vector<tuple<string, int>> found_files;
    string filename;
    string within;
    for (const auto & entry : fs::directory_iterator(lookup_dir)){
        filename = split(entry.path().string(), "/").back();
        if(fits_pattern(filename, pattern, &within)){
            found_files.push_back(tuple<string, int>{filename, stoi(within)});
            }
        }
    sort(found_files.begin(), found_files.end(), [=](tuple<string, int>& tup1, tuple<string, int>& tup2){return get<int>(tup1) < get<int>(tup2);});
    vector<string> sorted_files;
    for(tuple<string, int>& entry : found_files){
        sorted_files.push_back(get<string>(entry));
        }
    return sorted_files;
    }

array<float, 3> Helper::calc_cell_spans(array<float, 3> box, int num_splits){
    array<float, 3> result;
    float max_size = *max_element(box.begin(), box.end());
    float preferred_split = max_size/num_splits;
    int span;
    for(int i = 0; i < 3;  i++){
        span = round(box[i]/preferred_split);
        result[i] = box[i]/(span ? span : 1);
        }
    return result;
    }

int Helper::true_modulo(int divident, int divisor){
    if(divident < 0){
        return divisor + divident;
        }
    else{
        return divident%divisor;
        }
    }
