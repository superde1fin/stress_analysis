#pragma once

#include <array>

#include "Grid.hpp"

using namespace std;

class MaskGrid: public Grid{
    public:
        MaskGrid(array<float, 3> box, int num_splits);
        void add_mask(array<int, 3>);
        bool in_mask(array<float, 3> position);
    };
