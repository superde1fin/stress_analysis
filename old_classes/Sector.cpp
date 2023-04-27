#include <cmath>
#include <vector>
#include <array>
#include <iostream>

#include "Sector.hpp"
#include "Ellipse.hpp"
#include "Atom.hpp"
#include "Helper.hpp"

vector<array<float, 2>> Sector::segment_union;

Sector::Sector(Ellipse crcl, Atom atm, bool* success_indicator){
    Sector::thickness = crcl.get_thickness();
    float half_angle = Sector::get_half_angle(crcl, atm);
    float angular_position = Sector::get_angular_position(crcl, atm);
    lower_bound = angular_position - half_angle;
    upper_bound = angular_position + half_angle;
    /*
    cout << atm.get_id() << "---------------\n";
    cout << "Atom: " << atm.get_y() << " " << atm.get_z() << endl;
    cout << "Center: " << crcl.get_y() << " " << crcl.get_z() << endl;
    cout << "Bounds: " << lower_bound*180/M_PI << " " << upper_bound*180/M_PI<< endl;
    cout << "Axies: " << crcl.get_minor() << " " << crcl.get_major() << endl;
    */
    if(Sector::in_segments(crcl, atm)){
        *success_indicator = false;
        }else{
            *success_indicator = true;
            Sector::add2union();
            }
//    cout << "Area: " << Sector::get_area(crcl) << endl;
    }

bool Sector::in_segments(Ellipse crcl, Atom atm){
    float angle = Sector::get_angular_position(crcl, atm);
    int segments_num = Sector::segment_union.size();
    for(int i = 0; i < segments_num; i++){
        if(segment_union[i][0] <= angle && angle <= segment_union[i][1]){return true;}
        }
    return false;
    }

float Sector::get_half_angle(Ellipse elp, Atom atm){
    float R = pow(pow(atm.get_y() - elp.get_y(), 2) + pow((atm.get_z() - elp.get_z()), 2), 0.5); //Distance from the center of ellipse to the center of the atom
    return asin(atm.get_radius()/R);
    }
float Sector::get_angular_position(Ellipse elp, Atom atm){
    float theta = atan((atm.get_z() - elp.get_z())/(atm.get_y() - elp.get_y()));
    if(atm.get_y() == elp.get_y()){return M_PI/2*(atm.get_z() > elp.get_z() ? 1 : -1);}
    if(atm.get_y() < elp.get_y()){return theta + M_PI;}
    return theta;
    }
void Sector::add2union(){
    int segments_num = Sector::segment_union.size();
    bool intersection = false;
    if(!segments_num){
        Sector::segment_union.push_back(array<float, 2>{Sector::lower_bound, Sector::upper_bound});
        }else{
            for(int i = 0; i < segments_num; i++){
                if(Sector::lower_bound <= Sector::segment_union[i][0] && Sector::segment_union[i][0] <= Sector::upper_bound){
                    Sector::segment_union[i][0] = Sector::lower_bound;
                    intersection = true;
                    }
                if(Sector::lower_bound <= Sector::segment_union[i][1] && Sector::segment_union[i][1] <= Sector::upper_bound){
                    Sector::segment_union[i][1] = Sector::upper_bound;
                    intersection = true;
                    }
                }
            if(!intersection){
                Sector::segment_union.push_back(array<float, 2>{Sector::lower_bound, Sector::upper_bound});
                }
            }
    Helper::remove_dupl(segment_union);
    }

float Sector::get_area(Ellipse crcl){
    float major = crcl.get_major();
    float minor = crcl.get_minor();
    float area = 0;
    int segments_num = Sector::segment_union.size();
    for(int i = 0; i < segments_num; i++){
        float theta2 = segment_union[i][1];
        float theta1 = segment_union[i][0];
        area += 0.5*major*minor*(atan((major/minor)*tan(theta2)) - atan((major/minor)*tan(theta1)));
        //area += ((segment_union[i][1] - segment_union[i][0])/2)*pow(crcl.get_radius(), 2);
        }
    return area;
    }

void Sector::reset(){
    Sector::segment_union.clear();
    }

