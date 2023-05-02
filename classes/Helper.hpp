#pragma once

#include <string>
#include <vector>
#include <array>
#include <filesystem>
#include <tuple>

#include "Atom.hpp"

using namespace std;

class Helper{
    public:
        static void remove_dupl(std::vector<Atom> &v);
        static vector<string> split(string str, string delim);
        static bool string_contains(string to_check, string substr);
        static float dist(array<float, 3>atom1, array<float, 3> atom2, array<float, 3> box);
        static tuple<Atom, float> find_closest(Atom atm, vector<Atom> to_compare, array<float, 3> box);
        static void vector2d_csv(string name, string first_line, vector<vector<float>> vect);
        template <typename T>
        static string to_str(T value ){
            ostringstream ss;
            ss << value;
            return ss.str();
            }
        static vector<string> files_by_pattern(string cwd, string pattern);
        static vector<string> files_by_pattern(string cwd, string pattern, bool sort);
        static bool fits_pattern(string to_check, string pattern);
        static bool fits_pattern(string to_check, string pattern, string* filler);
        static tuple<array<float, 3>, array<int, 3>> calc_cell_spans(array<float, 3> box, int num_splits);
        static int true_modulo(int divident, int divisor);
    };
