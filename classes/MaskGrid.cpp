#include <array>

#include "MaskGrid.hpp"
#include "Grid.hpp"
#include "GridCell.hpp"

using namespace std;

MaskGrid::MaskGrid(array<float, 3> box, int num_splits) : Grid(box, num_splits) {}

// Adds a masked grid cell
void MaskGrid::add_mask(array<int, 3> key) {
    // If not masked already
    if (!MaskGrid::grid_map[key]) {
        // Create new cell
        GridCell* gc_ptr = new GridCell();
        // Add to mask
        MaskGrid::grid_map[key] = gc_ptr;
    }
}

// Checks whether a specific key is in the mask
bool MaskGrid::in_mask(array<float, 3> position) {
    array<int, 3> key = { static_cast<int>(position[0] / MaskGrid::float_sides[0]),
                          static_cast<int>(position[1] / MaskGrid::float_sides[1]),
                          static_cast<int>(position[2] / MaskGrid::float_sides[2]) };
    return MaskGrid::grid_map.count(key) > 0;
}
