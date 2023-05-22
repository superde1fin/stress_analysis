#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "Molecule.hpp"
#include "Helper.hpp"

using namespace std;

Molecule::Molecule(int type, string type_str, int id){
    Molecule::type = type;
    Molecule::type_str = type_str;
    Molecule::id = id;
    }

vector<Molecule*> Molecule::get_bonds(){
    return Molecule::bonds;
    }

Molecule* Molecule::add_bond(int type, string type_str, int id, int dist){
    Molecule* ml = new Molecule(type, type_str, id);
    Molecule::add_bond(ml, dist, 0);
    return ml;
    }

Molecule* Molecule::add_bond(Molecule* ml, int dist){
    return Molecule::add_bond(ml, dist, 0);
    }

Molecule* Molecule::add_bond(Molecule* ml, int dist, int rec_ctr){
    rec_ctr++;
    int num_bonds = Molecule::bonds.size(), i = 0;
    while(i < num_bonds && dist < Molecule::distances[i]){
        i++;
        }
    Molecule::bonds.insert(next(Molecule::bonds.begin(), i), ml);
    Molecule::distances.insert(next(Molecule::distances.begin(), i), dist);
    if(rec_ctr < 2){
        ml -> add_bond(this, dist, rec_ctr);
        }
    return ml;
    }

bool Molecule::is_atom(){
    return !Molecule::bonds.size();
    }

bool Molecule::in_molecule(int type){
    if(type == Molecule::type){
        return true;
        }
    else{
        bool res = false;
        for(Molecule* ml : Molecule::bonds){
            res = ml -> in_molecule(type);
            if(res){break;}
            }
        return res;
        }
    }

Molecule::~Molecule(){
    for(Molecule* ml : Molecule::bonds){
        delete ml;
        }
    }

map<int, tuple<vector<int>, int, string>> Molecule::get_count(set<int>* exclude){
//    cout << "Looking at molecule: " << Molecule::type << " with id: " << Molecule::id << endl;
    bool done;
    map<int, tuple<vector<int>, int, string>> count;
    vector<int> ids = {Molecule::id};
    get<vector<int>>(count[Molecule::type]) = ids;
    get<int>(count[Molecule::type]) = 1;
    get<string>(count[Molecule::type]) = Molecule::type_str;
    map<int, tuple<vector<int>, int, string>> temp_count;
    vector<int> non_intersection;
    vector<int>* existing_ids;
    exclude -> insert(Molecule::id);
    for(Molecule* ml : Molecule::bonds){
        if(!Helper::element_in(ml -> id, *exclude)){
            if(ml -> is_atom()){
//                cout << "The molecule " << ml -> id << " is an atom\n";
                get<vector<int>>(count[ml -> get_type()]).push_back(ml -> id);
                get<int>(count[ml -> get_type()])++;
                get<string>(count[ml -> get_type()]) = ml -> get_type_str();
//                cout << ml -> get_type() << " " << get<int>(count[ml -> get_type()]) << " " << get<string>(count[ml -> get_type()]) << endl;
                }
            else{
                temp_count = ml -> get_count(exclude);
                for(auto it = temp_count.begin(); it != temp_count.end(); ++it){
                    if(count.count(it -> first)){
                        non_intersection = Helper::not_in(get<vector<int>>(count[it -> first]), get<vector<int>>(it -> second));
                        get<int>(count[it -> first]) += non_intersection.size();
                        existing_ids = &get<vector<int>>(count[it -> first]);
                        (*existing_ids).insert(existing_ids -> end(), non_intersection.begin(), non_intersection.end());
                        get<string>(count[it -> first]) = get<string>(it -> second);
                        }
                    else{
                        get<int>(count[it -> first]) = get<int>(it -> second);
                        get<vector<int>>(count[it -> first]) = get<vector<int>>(it -> second);
                        get<string>(count[it -> first]) = get<string>(it -> second);
                        }
                    }
                }
            }
        }
//    cout << "----------------- Current count for atom: " << Molecule::id << endl;
//    for(auto it = count.begin(); it != count.end(); ++it){
//        cout << it -> first << " " << get<string>(it -> second) << " " << get<int>(it -> second) << endl;
//        }
    return count;
    }

string Molecule::get_name(){
    string name = "";
    set<int> exclude;
    map<int, tuple<vector<int>, int, string>> count = Molecule::get_count(&exclude);
    for(auto it = count.begin(); it != count.end(); ++it){
        name += (get<string>(it -> second) + "(" + Helper::to_str(get<int>(it -> second)) + ")");
        }
    return name;
    }


string Molecule::get_type_str(){
    return Molecule::type_str;
    }

int Molecule::get_type(){
    return Molecule::type;
    }

int Molecule::count_atoms(){
    int total = 0;
    set<int> exclude;
    map<int, tuple<vector<int>, int, string>> count = Molecule::get_count(&exclude);
    for(auto it = count.begin(); it != count.end(); ++it){
        total += get<int>(it -> second);
        }
    return total;
    }

tuple<int, Molecule*> Molecule::closest(int type, set<int> *exclude){
    tuple<int, Molecule*> res;
    get<int>(res) = numeric_limits<int>::max();
    get<Molecule*>(res) = this;
    tuple<int, Molecule*> iter_res;
    if(Molecule::get_type() == type){
        get<int>(res) = 0;
        return res;
        }
    else{
        if(Molecule::is_atom()){
            exclude -> insert(Molecule::id);
            return res;
            }
        }
    exclude -> insert(Molecule::id);
    for(Molecule* ml : Molecule::bonds){
        if(!Helper::element_in(ml -> id, *exclude)){
            iter_res = ml -> closest(type, exclude);
            if(get<int>(iter_res) < get<int>(res)){
                res = iter_res;
                }
            }
        }
    if(get<int>(res) != numeric_limits<int>::max()){
        get<int>(res)++;
        }
    return res;
    }

tuple<int, Molecule*> Molecule::closest(int type){
    set<int> exclude;
    return Molecule::closest(type, &exclude);
    }

Molecule* Molecule::truncate(int depth, Molecule* parent, int charge){
    charge += Molecule::charges[Molecule::type];
    Molecule* result = Molecule::copy();
    if(depth - 1){
        float dist;
        int bonds_size = Molecule::bonds.size();
        Molecule* ml;
        for(int i = 0; i < bonds_size; i++){
            ml = Molecule::bonds[i];
            dist = Molecule::distances[i];
            if(ml != parent){
                if(charge + Molecule::charges[ml -> get_type()] < 2){
                    result -> add_bond(ml -> truncate(depth - 1, this, charge), dist);
                    }
                else{
                    break;
                    }
                }
            }
        }
    return result;
    }

Molecule* Molecule::truncate(int depth){
    return Molecule::truncate(depth, this, 0);
    }

vector<int> Molecule::get_ids(){
    vector<int> result;
    vector<int> ids;
    set<int> exclude;
    map<int, tuple<vector<int>, int, string>> count = Molecule::get_count(&exclude);
    for(auto it = count.begin(); it != count.end(); ++it){
        ids = get<vector<int>>(it -> second);
        result.insert(result.end(), ids.begin(), ids.end());
        }
    return result;
    }

Molecule* Molecule::copy(){
    return new Molecule(Molecule::type, Molecule::type_str, Molecule::id);
    }

int Molecule::get_id(){
    return Molecule::id;
    }
