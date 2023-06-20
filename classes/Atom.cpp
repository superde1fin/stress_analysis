#include <iostream>
#include <cmath>
#include "Atom.hpp"

#include "Helper.hpp"

using namespace std;

//Setter functions
void Atom::set_id(int id){Atom::id = id;}
void Atom::set_pure_line(string line){Atom::pure_line = line;}
void Atom::set_type(int type){Atom::type = type;}
void Atom::set_position(array<float, 3> coordinates){
    Atom::x = coordinates[0];
    Atom::y = coordinates[1];
    Atom::z = coordinates[2];
    }
void Atom::set_stress_tensor(array<float, 6> stress_tensor){
    Atom::stress_tensor = stress_tensor;
    }
void Atom::set_radius(float radius){
    Atom::radius = radius;
    }
void Atom::set_potential(float potential){
    Atom::potential = potential;
    }
//Getter functions
int Atom::get_id(){return Atom::id;}
float Atom::get_potential(){return Atom::potential;}
int Atom::get_type(){return Atom::type;}
int Atom::get_radius(){return Atom::radius;}
array<float, 3> Atom::get_position(){return array<float, 3>{Atom::x, Atom::y, Atom::z};}
float Atom::get_x(){return Atom::x;}
float Atom::get_y(){return Atom::y;}
float Atom::get_z(){return Atom::z;}
string Atom::get_pure_line(){return Atom::pure_line;}
float Atom::get_ave_stress(){
    return pow(pow(Atom::stress_tensor[0], 2) + pow(Atom::stress_tensor[1], 2) + pow(Atom::stress_tensor[2], 2), 0.5);
    }
float Atom::get_stress_comp(int index){
    return Atom::stress_tensor[index];
    }
