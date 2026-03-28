#include "particles/container/cells/Cell.h"

#include <spdlog/spdlog.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace mol_sim;

void Cell::updateCache() {
    if (cache_dirty) {
        sorted_cache.clear();
        sorted_cache.reserve(indices.size());
        sorted_cache.insert(sorted_cache.end(), indices.begin(), indices.end());
        std::sort(sorted_cache.begin(), sorted_cache.end());  // NOLINT
        cache_dirty = false;
    }
}

Cell::Cell(CellType cell_type, std::array<double, 6> bounds) noexcept : type(cell_type), bounds(bounds) {}

void Cell::addParticle(size_t idx) {
    auto [_, inserted] = indices.insert(idx);
    cache_dirty = inserted || cache_dirty;
}

void Cell::removeParticle(size_t idx) noexcept { cache_dirty = indices.erase(idx) != 0 || cache_dirty; }

void Cell::updateParticleIndex(size_t old_idx, size_t new_idx) {
    auto it = indices.find(old_idx);
    if (it != indices.end()) {
        indices.erase(it);
        indices.insert(new_idx);
        cache_dirty = true;
    }
}

void Cell::clear() noexcept { indices.clear(); }
bool Cell::fits(R3 x) const noexcept {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t Cell::size() const noexcept { return indices.size(); }