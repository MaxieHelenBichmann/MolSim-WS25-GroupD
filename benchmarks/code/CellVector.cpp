#include "CellVector.h"

#include <algorithm>

using namespace mol_sim;

CellVector::CellVector(std::array<double, 6> bounds) : bounds(bounds) {}
void CellVector::addParticle(size_t idx) {
    auto it = std::find(indices.begin(), indices.end(), idx);  // NOLINT
    if (it == indices.end()) {
        indices.push_back(idx);
    }
}
void CellVector::removeParticle(size_t idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == idx) {
            indices.erase(it);
            return;
        }
    }
}
void CellVector::updateParticleIndex(size_t old_idx, size_t new_idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == old_idx) {
            indices.erase(it);
            indices.push_back(new_idx);
            return;
        }
    }
}
void CellVector::clear() { indices.clear(); }
std::vector<size_t>& CellVector::particles() { return indices; }
const std::vector<size_t>& CellVector::particles() const { return indices; }
bool CellVector::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t CellVector::size() { return indices.size(); }