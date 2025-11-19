#include "CellSet.h"

using namespace mol_sim;

CellSet::CellSet(std::array<double, 6> bounds) : bounds(bounds) {}
void CellSet::addParticle(size_t idx) { indices.insert(idx); }
void CellSet::removeParticle(size_t idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == idx) {
            indices.erase(it);
            return;
        }
    }
}
void CellSet::updateParticleIndex(size_t old_idx, size_t new_idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == old_idx) {
            indices.erase(it);
            indices.insert(new_idx);
            return;
        }
    }
}
void CellSet::clear() { indices.clear(); }
std::set<size_t>& CellSet::particles() { return indices; }
const std::set<size_t>& CellSet::particles() const { return indices; }
bool CellSet::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t CellSet::size() { return indices.size(); }