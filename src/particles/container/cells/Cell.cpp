#include "particles/container/cells/Cell.h"

#include <spdlog/spdlog.h>

#include <algorithm>

using namespace mol_sim;

Cell::Cell(CellType cell_type, std::array<double, 6> bounds) : type(cell_type), bounds(bounds) {}

void Cell::addParticle(size_t idx) { indices.insert(idx); }
void Cell::removeParticle(size_t idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == idx) {
            indices.erase(it);  // NOLINT
            return;
        }
    }
}
void Cell::updateParticleIndex(size_t old_idx, size_t new_idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == old_idx) {
            indices.erase(it);
            indices.insert(new_idx);
            return;
        }
    }
}
void Cell::clear() { indices.clear(); }
std::set<size_t>& Cell::particles() { return indices; }
const std::set<size_t>& Cell::particles() const { return indices; }
bool Cell::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t Cell::size() { return indices.size(); }
CellType Cell::getType() { return type; }