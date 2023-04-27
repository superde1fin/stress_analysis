#include <array>
#include <limits>
#include <vector>

#include "Bubble.hpp"
#include "Atom.hpp"
#include "Helper.hpp"

using namespace std;

Bubble::Bubble(array<float, 3> position, array<float, 3> span){
    Bubble::x = position[0];
    Bubble::y = position[1];
    Bubble::z = position[2];
    Bubble::x_span = span[0];
    Bubble::y_span = span[1];
    Bubble::z_span = span[2];
    Bubble::position = position;
    Bubble::span = span;
    }

array<int, 3> Bubble::get_mapping(){
    return array<int, 3>{(int)(Bubble::x/Bubble::x_span), (int)(Bubble::y/Bubble::y_span), (int)(Bubble::z/Bubble::z_span)};
    }

array<float, 3> Bubble::get_position(){
    return Bubble::position;
    }

array<float, 3> Bubble::get_span(){
    return Bubble::span;
    }
