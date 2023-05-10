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
#include "Molecule.hpp"

using namespace std;
namespace fs = filesystem;

//Removes the duplicates from a vector of atoms by ids
void Helper::remove_dupl(vector<Atom>& atoms){
    map<int, int> atom_ctrs;
    vector<Atom> result;
    //Loop through all atoms
    for(auto it = atoms.begin(); it != atoms.end(); it++){
        //If has not been recorded before
        if(!atom_ctrs[(*it).get_id()]){
            //Add to a final vector
            result.push_back(*it);
            }
        //Record an occurance of the atom
        atom_ctrs[(*it).get_id()]++;
        }
    //Reassign the atoms vector by reference
    atoms = result;
    }

//Split a line by delimeter
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

//Check whether a string contains a substring
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

//Calculate the distance between two points in a periodic box
float Helper::dist(array<float, 3>atom1, array<float, 3> atom2, array<float, 3> box){
    float x, y, z;
    x = abs(atom1[0] - atom2[0]);
    x = (x > box[0]/2) ? (box[0] - x) : x;
    
    y = abs(atom1[1] - atom2[1]);
    y = (y > box[1]/2) ? (box[1] - y) : y;
    
    z = abs(atom1[2] - atom2[2]);
    z = (z > box[2]/2) ? (box[2] - z) : z;
    
    return pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5);
    }

//Returns closest atom out of th given vector
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

//Record the contents of a 2d vector in a csv file
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


//Check whether a string fits a specific pattern
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

//The same function as above but returns what was found in place of *
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

//Returns a vector of filenames that fit the pattern
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

//The same function as above, but sorts the results by filler
vector<string> Helper::files_by_pattern(string lookup_dir, string pattern, vector<string>* pattern_fillers, bool sort_results = false){
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
        pattern_fillers -> push_back(Helper::to_str(get<int>(entry)));
        }
    return sorted_files;
    }

//Calculates appropriate grid size based on the number of splits provided
tuple<array<float, 3>, array<int, 3>> Helper::calc_cell_spans(array<float, 3> box, int num_splits){
    array<float, 3> float_sides;
    array<int, 3> int_sides;
    float max_size = *max_element(box.begin(), box.end());
    float preferred_split = max_size/num_splits;
    int span;
    for(int i = 0; i < 3;  i++){
        span = round(box[i]/preferred_split);
        span = span ? span : 1;
        int_sides[i] = span;
        float_sides[i] = box[i]/(span);
        }
    return make_tuple(float_sides, int_sides);
    }

//Implements the correct modulo that works with negative numbers (small)
int Helper::true_modulo(int divident, int divisor){
    if(divident < 0){
        return divisor + divident;
        }
    else{
        return divident%divisor;
        }
    }

bool Helper::element_in(int elem, set<int> search_vector){
    for(int item : search_vector){
        if(item == elem){
            return true;
            }
        }
    return false;
    }

void Helper::vector_of_maps2csv(string name, vector<map<string, float>> map_vect, vector<string> first_col){
    cout << "Writing to CSV\n";
    map<string, int> first_line;
    int number_cols = 1;
    first_line["Timestep"] = 0;
    
    //Record all the columns that are going to be present in a csv
    for(map<string, float> line: map_vect){
        for(map<string, float>::iterator it_sp = line.begin(); it_sp != line.end(); ++it_sp){
            if(!first_line.count(it_sp -> first)){
                first_line[it_sp -> first] =  number_cols++;
                }
            }
        }
        
    //Create a csv body string
    int size_line = first_line.size();
    vector<string> line_vect(size_line, "0");
    int vsize = map_vect.size();
    int fsize = first_col.size();
    int size = (vsize > fsize) ? fsize : vsize;
    map<string, float> line;
    string csv_body("");
    
    for(int i = 0; i < size; i++){
        line = map_vect[i];
        line_vect[0] = first_col[i];
        for(map<string, float>::iterator it_sp = line.begin(); it_sp != line.end(); ++it_sp){
            line_vect[first_line[it_sp -> first]] = to_str(it_sp -> second);
            }
        for(int j = 0; j < size_line; j++){
            csv_body += to_str(line_vect[j]) + ((j == size_line - 1) ? "\n" : ", ");
            }
        }
        
    string first_line_str("");
    map<int, string> reversed_fl;
    for(map<string, int>::iterator it_name = first_line.begin(); it_name != first_line.end(); ++it_name){
        reversed_fl[it_name -> second] = it_name -> first;
        }
        
        
    for(int i = 0; i < size_line; i++){
        first_line_str += reversed_fl[i] + ((i == size_line - 1) ? "\n" : ", ");
        }
    
    
    ofstream myfile;
    myfile.open(name + ".csv");
    myfile << first_line_str << csv_body;
    myfile.close();
   myfile.close();
    }

bool Helper::element_in(int elem, vector<int> search_vector){
    for(int item : search_vector){
        if(item == elem){
            return true;
            }
        }
    return false;
    }

vector<int> Helper::not_in(vector<int> a, vector<int> b){
    vector<int> non_intersection;
    for(int item : b){
        if(!Helper::element_in(item, a)){
            non_intersection.push_back(item);
            }
        }
    return non_intersection;
    }
