#include "particles/container/cells/Cell.h"

#include <spdlog/spdlog.h>

#include <algorithm>

using namespace mol_sim;

Cell::Cell(CellType cell_type, std::array<double, 6> bounds) : type(cell_type), bounds(bounds) {}

void Cell::addParticle(size_t idx) {
    if (std::find(indices.begin(), indices.end(), idx) == indices.end()) {  // NOLINT
        indices.push_back(idx);
    }
}
void Cell::removeParticle(size_t idx) {
    for (size_t i = 0; i < indices.size(); ++i) {  // NOLINT
        if (indices[i] == idx) {
            indices.erase(indices.begin() + i);  // NOLINT
            return;
        }
    }
}
void Cell::updateParticleIndex(size_t old_idx, size_t new_idx) {
    for (size_t i = 0; i < indices.size(); ++i) {  // NOLINT
        if (indices[i] == old_idx) {
            indices[i] = new_idx;
            return;
        }
    }
}
void Cell::clear() { indices.clear(); }
std::vector<size_t>& Cell::particles() { return indices; }
const std::vector<size_t>& Cell::particles() const { return indices; }
size_t Cell::operator[](size_t idx) { return indices[idx]; }
bool Cell::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t Cell::size() { return indices.size(); }
CellType Cell::getType() { return type; }