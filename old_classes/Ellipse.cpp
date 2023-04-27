#include <cmath>
#include <array>
#include <iostream>

#include "Ellipse.hpp"
#include "Atom.hpp"
#include "Helper.hpp"

Ellipse::Ellipse(array<float, 3> origin, float thickness, Atom atm){
    Ellipse::x = origin[0];
    Ellipse::y = origin[1];
    Ellipse::z = origin[2];
    float axes_ratio = 5.0;
    Ellipse::minor = pow(pow(atm.get_y() - Ellipse::y, 2) + pow(axes_ratio, 2)*pow(atm.get_z() - Ellipse::z, 2), 0.5)/axes_ratio;
    Ellipse::major = 2*minor;
    }
bool Ellipse::within(Atom atm){
    return pow(atm.get_y() - Ellipse::y, 2)/pow(Ellipse::major, 2) + pow(atm.get_z() - Ellipse::z, 2)/pow(Ellipse::minor, 2) <= 1 && Ellipse::x - (Ellipse::thickness)/2 < atm.get_x() < Ellipse::x + (Ellipse::thickness)/2;
    }
float Ellipse::get_area(){
    return M_PI*Ellipse::major*Ellipse::minor;
    }
float Ellipse::get_x(){
    return Ellipse::x;
    }
float Ellipse::get_y(){
    return Ellipse::y;
    }
float Ellipse::get_z(){
    return Ellipse::z;
    }
float Ellipse::get_minor(){
    return Ellipse::minor;
    }
float Ellipse::get_major(){
    return Ellipse::major;
    }
float Ellipse::get_thickness(){
    return Ellipse::thickness;
    }
