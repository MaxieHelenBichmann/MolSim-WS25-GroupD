#include "CellUnorderedSet.h"

using namespace mol_sim;

CellUnorderedSet::CellUnorderedSet(std::array<double, 6> bounds) : bounds(bounds) {}
void CellUnorderedSet::addParticle(size_t idx) { indices.insert(idx); }
void CellUnorderedSet::removeParticle(size_t idx) {
    auto it = indices.find(idx);
    if (it != indices.end()) {
        indices.erase(it);
    }
}
void CellUnorderedSet::updateParticleIndex(size_t old_idx, size_t new_idx) {
    auto it = indices.find(old_idx);
    if (it != indices.end()) {
        indices.erase(it);
        indices.insert(new_idx);
    }
}
void CellUnorderedSet::clear() { indices.clear(); }
std::unordered_set<size_t>& CellUnorderedSet::particles() { return indices; }
const std::unordered_set<size_t>& CellUnorderedSet::particles() const { return indices; }
bool CellUnorderedSet::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t CellUnorderedSet::size() { return indices.size(); }