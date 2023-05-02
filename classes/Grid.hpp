#pragma once

#include <array>

using namespace std;

class Grid{
    public:
        Grid(array<float, 3> box, int num_splits);
        ~Grid();
        array<int, 3> get_int_sides();
        array<float, 3> get_float_sides();
        int get_size();
    protected:
        array<float, 3> box;
        array<float, 3> float_sides;
        array<int, 3> int_sides;
        struct ArrayHasher {
            size_t operator()(const array<int, 3>& a) const {
                size_t h = 0;

                for (auto e : a) {
                    h ^= hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2); 
                    }
                return h;
                }   
            };
        unordered_map<array<int, 3>, GridCell*, ArrayHasher> grid_map;
    };
