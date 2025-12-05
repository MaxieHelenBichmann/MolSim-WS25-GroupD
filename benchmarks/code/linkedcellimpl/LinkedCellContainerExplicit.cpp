#include "LinkedCellContainerExplicit.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

using namespace mol_sim;

// --------------------------------------------------------------------------------------------------------------

CellExplicit::CellExplicit(CellType cell_type, std::array<double, 6> bounds) : bounds(bounds), type(cell_type) {}

void CellExplicit::addParticle(size_t idx) {
    auto it = std::find(indices.begin(), indices.end(), idx);  // NOLINT
    if (it == indices.end()) {
        indices.push_back(idx);
    }
}
void CellExplicit::removeParticle(size_t idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == idx) {
            indices.erase(it);
            return;
        }
    }
}
void CellExplicit::updateParticleIndex(size_t old_idx, size_t new_idx) {
    for (auto it = indices.begin(); it != indices.end(); ++it) {  // NOLINT
        if (*it == old_idx) {
            indices.erase(it);
            indices.push_back(new_idx);
            return;
        }
    }
}
void CellExplicit::clear() { indices.clear(); }
std::vector<size_t>& CellExplicit::particles() { return indices; }
const std::vector<size_t>& CellExplicit::particles() const { return indices; }
bool CellExplicit::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t CellExplicit::size() { return indices.size(); }
CellType CellExplicit::getType() { return type; }

// --------------------------------------------------------------------------------------------------------------

LinkedCellContainerExplicit::LinkedCellContainerExplicit(R3 domain_size, double cutoff_radius)
    : domain_size(domain_size), cutoff_radius(cutoff_radius) {
    for (size_t dim = 0; dim < 3; ++dim) {
        size_t inner_cells = 0U;
        if (cutoff_radius > 0.0) {
            inner_cells = static_cast<size_t>(std::floor(domain_size[dim] / cutoff_radius));
        }
        if (inner_cells == 0U) {
            inner_cells = 1U;
        }
        cell_length[dim] = domain_size[dim] / static_cast<double>(inner_cells);
        num_cells[dim] = inner_cells + 2U;
    }

    SPDLOG_DEBUG("Initializing Linked-Cell Container with dimensions: {} x {} x {}", num_cells[0], num_cells[1],
                 num_cells[2]);
    SPDLOG_DEBUG("and cell lengths: {} x {} x {}", cell_length[0], cell_length[1], cell_length[2]);

    cells.reserve(num_cells[0] * num_cells[1] * num_cells[2]);
    for (size_t z = 0; z < num_cells[2]; ++z) {
        for (size_t y = 0; y < num_cells[1]; ++y) {
            for (size_t x = 0; x < num_cells[0]; ++x) {
                const double x_min = (static_cast<double>(x) - 1.0) * cell_length[0];
                const double y_min = (static_cast<double>(y) - 1.0) * cell_length[1];
                const double z_min = (static_cast<double>(z) - 1.0) * cell_length[2];

                std::array<double, 6> bounds = {x_min, x_min + cell_length[0], y_min, y_min + cell_length[1],
                                                z_min, z_min + cell_length[2]};

                CellType type = CellType::INNER;
                const bool halo = (x == 0U || x == num_cells[0] - 1U || y == 0U || y == num_cells[1] - 1U || z == 0U ||
                                   z == num_cells[2] - 1U);
                if (halo) {
                    type = CellType::HALO;
                } else {
                    const bool boundary = (x == 1U || x == num_cells[0] - 2U || y == 1U || y == num_cells[1] - 2U ||
                                           z == 1U || z == num_cells[2] - 2U);
                    if (boundary) {
                        type = CellType::BOUNDARY;
                    }
                }

                cells.emplace_back(type, bounds);
            }
        }
    }
}

size_t LinkedCellContainerExplicit::findCellIndex(R3 vec) const {
    const auto index_for_dim = [this](double coord, size_t dim) {
        const double normalized = (coord + cell_length[dim]) / cell_length[dim];

        if (normalized <= 0.0) {
            return size_t{0};
        }

        if (normalized >= static_cast<double>(num_cells[dim] - 1U)) {
            return num_cells[dim] - 1U;
        }

        return static_cast<size_t>(normalized);
    };

    const size_t x_idx = index_for_dim(vec[0], 0U);
    const size_t y_idx = index_for_dim(vec[1], 1U);
    const size_t z_idx = index_for_dim(vec[2], 2U);

    return (z_idx * num_cells[1] * num_cells[0]) + (y_idx * num_cells[0]) + x_idx;
}

std::vector<CellExplicit*> LinkedCellContainerExplicit::findAdjacentCells(size_t cell_idx) {
    std::vector<CellExplicit*> adjacent_cells;
    size_t z_idx = cell_idx / (num_cells[0] * num_cells[1]);
    size_t y_idx = (cell_idx / num_cells[0]) % num_cells[1];
    size_t x_idx = cell_idx % num_cells[0];

    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                size_t nx = x_idx + dx;
                size_t ny = y_idx + dy;
                size_t nz = z_idx + dz;

                if (nx < num_cells[0] && ny < num_cells[1] && nz < num_cells[2]) {
                    size_t neighbor_idx = (nz * num_cells[1] * num_cells[0]) + (ny * num_cells[0]) + nx;
                    adjacent_cells.push_back(&cells[neighbor_idx]);
                }
            }
        }
    }

    return adjacent_cells;
}
std::vector<const CellExplicit*> LinkedCellContainerExplicit::findAdjacentCells(size_t cell_idx) const {
    std::vector<const CellExplicit*> adjacent_cells;
    size_t z_idx = cell_idx / (num_cells[0] * num_cells[1]);
    size_t y_idx = (cell_idx / num_cells[0]) % num_cells[1];
    size_t x_idx = cell_idx % num_cells[0];

    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                size_t nx = x_idx + dx;
                size_t ny = y_idx + dy;
                size_t nz = z_idx + dz;

                if (nx < num_cells[0] && ny < num_cells[1] && nz < num_cells[2]) {
                    size_t neighbor_idx = (nz * num_cells[1] * num_cells[0]) + (ny * num_cells[0]) + nx;
                    adjacent_cells.push_back(&cells[neighbor_idx]);
                }
            }
        }
    }

    return adjacent_cells;
}

void LinkedCellContainerExplicit::findBoundaryCells(const BoundaryLocation type,  // NOLINT
                                                    std::vector<const CellExplicit*>& boundary_cells,
                                                    size_t offset) const {
    switch (type) {
        case BoundaryLocation::UPPER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx =
                        ((num_cells[2] - (1U + offset)) * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::LOWER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (offset * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::LEFT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + offset;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::RIGHT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx =
                        (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + (num_cells[0] - (1U + offset));
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::FRONT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (offset * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::BACK: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx =
                        (z * num_cells[1] * num_cells[0]) + ((num_cells[1] - (1U + offset)) * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
    }
}
void LinkedCellContainerExplicit::findBoundaryCells(const BoundaryLocation type,  // NOLINT
                                                    std::vector<CellExplicit*>& boundary_cells, size_t offset) {
    switch (type) {
        case BoundaryLocation::UPPER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx =
                        ((num_cells[2] - (1U + offset)) * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::LOWER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (offset * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::LEFT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + offset;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::RIGHT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx =
                        (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + (num_cells[0] - (1U + offset));
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::FRONT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (offset * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryLocation::BACK: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx =
                        (z * num_cells[1] * num_cells[0]) + ((num_cells[1] - (1U + offset)) * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
    }
}
void LinkedCellContainerExplicit::decreaseCellIndices(size_t starting_idx) {
    for (size_t i = starting_idx; i < data.size(); ++i) {
        size_t cell_idx = findCellIndex(data[i].getX());
        cells[cell_idx].updateParticleIndex(i + 1, i);
    }
}

bool LinkedCellContainerExplicit::fitsDomain(R3 v) const {
    return (v[0] >= 0.0 && v[0] <= domain_size[0]) && (v[1] >= 0.0 && v[1] <= domain_size[1]) &&
           (v[2] >= 0.0 && v[2] <= domain_size[2]);
}

bool LinkedCellContainerExplicit::fitsContainer(R3 v) const {
    return (v[0] >= -cell_length[0] && v[0] <= domain_size[0] + cell_length[0]) &&
           (v[1] >= -cell_length[1] && v[1] <= domain_size[1] + cell_length[1]) &&
           (v[2] >= -cell_length[2] && v[2] <= domain_size[2] + cell_length[2]);
}

Particle& LinkedCellContainerExplicit::operator[](size_t idx) { return data[idx]; }
const Particle& LinkedCellContainerExplicit::operator[](size_t idx) const { return data[idx]; }

size_t LinkedCellContainerExplicit::size() const { return data.size(); }
bool LinkedCellContainerExplicit::empty() const { return data.empty(); }

// modify
void LinkedCellContainerExplicit::clear() {
    data.clear();
    for (auto& cell : cells) {
        cell.clear();
    }
}
void LinkedCellContainerExplicit::reserve(size_t n) { data.reserve(n); }

void LinkedCellContainerExplicit::addParticle(Particle&& value) {
    if (!fitsContainer(value.getX())) {
        return;
    }
    data.push_back(value);
    cells[findCellIndex(value.getX())].addParticle(data.size() - 1);
}
void LinkedCellContainerExplicit::addParticle(const Particle& value) {
    if (!fitsContainer(value.getX())) {
        return;
    }
    data.push_back(value);
    cells[findCellIndex(value.getX())].addParticle(data.size() - 1);
}

void LinkedCellContainerExplicit::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    if (!fitsContainer(x_arg)) {
        return;
    }
    data.emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
    cells[findCellIndex(x_arg)].addParticle(data.size() - 1);
}
void LinkedCellContainerExplicit::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg,
                                              int type) {
    if (!fitsContainer(x_arg)) {
        return;
    }
    data.emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
    cells[findCellIndex(x_arg)].addParticle(data.size() - 1);
}

std::vector<Particle>::iterator LinkedCellContainerExplicit::eraseParticle(std::vector<Particle>::iterator p) {
    size_t cell_idx = findCellIndex(p->getX());

    if (cell_idx < cells.size()) {
        size_t idx = &(*p) - data.data();
        cells[cell_idx].removeParticle(idx);
        std::vector<Particle>::iterator it = data.erase(p);  // NOLINT
        if (it != data.end()) {
            decreaseCellIndices(idx);
        }
        return it;
    }
    return data.end();
}

std::vector<Particle>::iterator LinkedCellContainerExplicit::updateParticlePosition(std::vector<Particle>::iterator p,
                                                                                    R3 new_x) {
    size_t old_cell_idx = findCellIndex(p->getX());

    if (!fitsContainer(new_x)) {
        // Particle moved completely outside container
        return eraseParticle(p);
    }

    size_t new_cell_idx = findCellIndex(new_x);

    if (new_cell_idx != old_cell_idx) {
        cells[old_cell_idx].removeParticle(static_cast<size_t>(p - data.begin()));
        cells[new_cell_idx].addParticle(static_cast<size_t>(p - data.begin()));
    }

    p->getX() = new_x;
    return ++p;
}

// normal iterators
std::vector<Particle>::iterator LinkedCellContainerExplicit::begin() { return data.begin(); }
std::vector<Particle>::const_iterator LinkedCellContainerExplicit::begin() const { return data.begin(); }
std::vector<Particle>::const_iterator LinkedCellContainerExplicit::cbegin() const { return data.cbegin(); }
std::vector<Particle>::iterator LinkedCellContainerExplicit::end() { return data.end(); }
std::vector<Particle>::const_iterator LinkedCellContainerExplicit::end() const { return data.end(); }
std::vector<Particle>::const_iterator LinkedCellContainerExplicit::cend() const { return data.cend(); }

// proximity iterators
LinkedCellContainerExplicit::proximity_iterator LinkedCellContainerExplicit::proximityBegin(
    R3 center, [[maybe_unused]] size_t offset) {
    std::vector<CellExplicit*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<CellExplicit*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (CellExplicit* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    if (nonempty_adjacent_cells.empty()) {
        std::vector<CellExplicit*> same_cell = {&cells[findCellIndex(center)]};
        return proximity_iterator{center, cutoff_radius, same_cell.front()->particles().begin(), same_cell, &data};
    }

    return proximity_iterator{center, cutoff_radius, nonempty_adjacent_cells.front()->particles().begin(),
                              nonempty_adjacent_cells, &data};
}

LinkedCellContainerExplicit::proximity_iterator LinkedCellContainerExplicit::proximityEnd(R3 center) {
    std::vector<CellExplicit*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<CellExplicit*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (CellExplicit* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    if (nonempty_adjacent_cells.empty()) {
        return proximityBegin(center);
    }

    return proximity_iterator{center, cutoff_radius, nonempty_adjacent_cells.back()->particles().end(),
                              nonempty_adjacent_cells, &data};
}

LinkedCellContainerExplicit::const_proximity_iterator LinkedCellContainerExplicit::proximityBegin(
    R3 center, [[maybe_unused]] size_t offset) const {
    std::vector<const CellExplicit*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<const CellExplicit*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (const CellExplicit* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    if (nonempty_adjacent_cells.empty()) {
        std::vector<const CellExplicit*> same_cell = {&cells[findCellIndex(center)]};
        return const_proximity_iterator{center, cutoff_radius, same_cell.front()->particles().begin(), same_cell,
                                        &data};
    }

    return const_proximity_iterator{center, cutoff_radius, nonempty_adjacent_cells.front()->particles().begin(),
                                    nonempty_adjacent_cells, &data};
}
LinkedCellContainerExplicit::const_proximity_iterator LinkedCellContainerExplicit::proximityEnd(R3 center) const {
    std::vector<const CellExplicit*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<const CellExplicit*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (const CellExplicit* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    if (nonempty_adjacent_cells.empty()) {
        return proximityBegin(center);
    }

    return const_proximity_iterator{center, cutoff_radius, nonempty_adjacent_cells.back()->particles().end(),
                                    nonempty_adjacent_cells, &data};
}

// boundary and halo iterators
LinkedCellContainerExplicit::proximity_iterator LinkedCellContainerExplicit::haloBegin(
    const std::set<BoundaryLocation>& boundary_types) {
    std::vector<CellExplicit*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::set<CellExplicit*> seen;
    std::vector<CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    if (unique_and_nonempty_cells.empty()) {
        return proximity_iterator{
            R3{}, std::numeric_limits<double>::infinity(), cells.front().particles().begin(), {&cells.front()}, &data};
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.front()->particles().begin(), unique_and_nonempty_cells, &data};
}
LinkedCellContainerExplicit::const_proximity_iterator LinkedCellContainerExplicit::haloBegin(
    const std::set<BoundaryLocation>& boundary_types) const {
    std::vector<const CellExplicit*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::set<const CellExplicit*> seen;
    std::vector<const CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    if (unique_and_nonempty_cells.empty()) {
        return const_proximity_iterator{
            R3{}, std::numeric_limits<double>::infinity(), cells.front().particles().begin(), {&cells.front()}, &data};
    }

    return const_proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                                    unique_and_nonempty_cells.front()->particles().begin(), unique_and_nonempty_cells,
                                    &data};
}
LinkedCellContainerExplicit::proximity_iterator LinkedCellContainerExplicit::haloEnd(
    const std::set<BoundaryLocation>& boundary_types) {
    std::vector<CellExplicit*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::unordered_set<CellExplicit*> seen;
    std::vector<CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    if (unique_and_nonempty_cells.empty()) {
        return haloBegin(boundary_types);
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.back()->particles().end(), unique_and_nonempty_cells, &data};
}
LinkedCellContainerExplicit::const_proximity_iterator LinkedCellContainerExplicit::haloEnd(
    const std::set<BoundaryLocation>& boundary_types) const {
    std::vector<const CellExplicit*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::unordered_set<const CellExplicit*> seen;
    std::vector<const CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    if (unique_and_nonempty_cells.empty()) {
        return haloBegin(boundary_types);
    }

    return const_proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                                    unique_and_nonempty_cells.back()->particles().end(), unique_and_nonempty_cells,
                                    &data};
}

LinkedCellContainerExplicit::proximity_iterator LinkedCellContainerExplicit::boundaryBegin(
    const std::set<BoundaryLocation>& boundary_types) {
    std::vector<CellExplicit*> relevant_cells;
    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<CellExplicit*> seen;
    std::vector<CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    if (unique_and_nonempty_cells.empty()) {
        return proximity_iterator{
            R3{}, std::numeric_limits<double>::infinity(), cells.front().particles().begin(), {&cells.front()}, &data};
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.front()->particles().begin(), unique_and_nonempty_cells, &data};
}
LinkedCellContainerExplicit::const_proximity_iterator LinkedCellContainerExplicit::boundaryBegin(
    const std::set<BoundaryLocation>& boundary_types) const {
    std::vector<const CellExplicit*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<const CellExplicit*> seen;
    std::vector<const CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    if (unique_and_nonempty_cells.empty()) {
        return const_proximity_iterator{
            R3{}, std::numeric_limits<double>::infinity(), cells.front().particles().begin(), {&cells.front()}, &data};
    }

    return const_proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                                    unique_and_nonempty_cells.front()->particles().begin(), unique_and_nonempty_cells,
                                    &data};
}
LinkedCellContainerExplicit::proximity_iterator LinkedCellContainerExplicit::boundaryEnd(
    const std::set<BoundaryLocation>& boundary_types) {
    std::vector<CellExplicit*> relevant_cells;
    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<CellExplicit*> seen;
    std::vector<CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }
    if (unique_and_nonempty_cells.empty()) {
        return boundaryBegin(boundary_types);
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.back()->particles().end(), unique_and_nonempty_cells, &data};
}
LinkedCellContainerExplicit::const_proximity_iterator LinkedCellContainerExplicit::boundaryEnd(
    const std::set<BoundaryLocation>& boundary_types) const {
    std::vector<const CellExplicit*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<const CellExplicit*> seen;
    std::vector<const CellExplicit*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellExplicit* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    if (unique_and_nonempty_cells.empty()) {
        return boundaryBegin(boundary_types);
    }

    return const_proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                                    unique_and_nonempty_cells.back()->particles().end(), unique_and_nonempty_cells,
                                    &data};
}

R3 LinkedCellContainerExplicit::getDomainSize() { return domain_size; }