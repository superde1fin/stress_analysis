#pragma once
#include <array>
#include <string>

using namespace std;

class Atom{
    private:
        int id;
        int type;
        float x;
        float y;
        float z;
        float potential;
        array<float, 6> stress_tensor;
        float radius;
        string pure_line;
        int Oid = 0;

        float get_ave_stress();

        //Setters
        void set_id(int id);
        void set_type(int type);
        void set_position(array<float, 3> coordinates);
        void set_stress_tensor(array<float, 6> stresss_tensor);
        void set_radius(float radius);
        void set_pure_line(string line);
        void set_potential(float potential);

    public:
        //Getters
        int get_id();
        int get_type();
        int get_radius();
        array<float, 3> get_position();
        float get_x();
        float get_y();
        float get_z();
        float get_stress_comp(int index);
        float get_potential();
        string get_pure_line();
    friend class System;
    friend class AtomGrid;
    };
