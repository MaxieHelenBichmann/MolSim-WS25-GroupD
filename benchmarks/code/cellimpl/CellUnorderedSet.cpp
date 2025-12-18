#include "CellUnorderedSet.h"

#include <algorithm>

using namespace mol_sim;

void CellUnorderedSet::updateCache() {
    if (cache_dirty) {
        sorted_cache.clear();
        sorted_cache.reserve(indices.size());
        sorted_cache.insert(sorted_cache.end(), indices.begin(), indices.end());
        std::sort(sorted_cache.begin(), sorted_cache.end());  // NOLINT
        cache_dirty = false;
    }
}

CellUnorderedSet::CellUnorderedSet(std::array<double, 6> bounds) : bounds(bounds) {}
void CellUnorderedSet::addParticle(size_t idx) {
    if (indices.insert(idx).second) {
        cache_dirty = true;
    }
}
void CellUnorderedSet::removeParticle(size_t idx) {
    auto it = indices.find(idx);
    if (it != indices.end()) {
        indices.erase(it);
        cache_dirty = true;
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

std::vector<size_t>::iterator CellUnorderedSet::stableIteratorBegin() {
    updateCache();
    return sorted_cache.begin();
}
std::vector<size_t>::iterator CellUnorderedSet::stableIteratorEnd() {
    updateCache();
    return sorted_cache.end();
}

bool CellUnorderedSet::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t CellUnorderedSet::size() { return indices.size(); }
