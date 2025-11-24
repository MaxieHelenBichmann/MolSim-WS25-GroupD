#include "LinkedCellContainerDirect.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

using namespace mol_sim;

// --------------------------------------------------------------------------------------------------------------

CellDirect::CellDirect(CellType cell_type, std::array<double, 6> bounds) : bounds(bounds), type(cell_type) {}

void CellDirect::addParticle(Particle&& value) {
    if (std::find(data.begin(), data.end(), value) == data.end()) {  // NOLINT
        data.push_back(value);
    }
}
void CellDirect::addParticle(const Particle& value) {
    if (std::find(data.begin(), data.end(), value) == data.end()) {  // NOLINT
        data.push_back(value);
    }
}

void CellDirect::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    if (std::find(data.begin(), data.end(), Particle(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg)) ==  // NOLINT
        data.end()) {
        data.emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
    }
}
void CellDirect::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type) {
    if (std::find(data.begin(), data.end(), Particle(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type)) ==  // NOLINT
        data.end()) {
        data.emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
    }
}
Particle CellDirect::removeParticle(size_t idx) {
    Particle p = data[idx];
    data.erase(data.begin() + idx);  // NOLINT
    return p;
}
void CellDirect::clear() { data.clear(); }
std::vector<Particle>& CellDirect::particles() { return data; }
const std::vector<Particle>& CellDirect::particles() const { return data; }
bool CellDirect::fits(R3 x) const {
    return bounds[0] <= x[0] && x[0] <= bounds[1] && bounds[2] <= x[1] && x[1] <= bounds[3] && bounds[4] <= x[2] &&
           x[2] <= bounds[5];
}

size_t CellDirect::size() { return data.size(); }
CellType CellDirect::getType() { return type; }

// --------------------------------------------------------------------------------------------------------------

LinkedCellContainerDirect::LinkedCellContainerDirect(R3 domain_size, double cutoff_radius) : domain_size(domain_size) {
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

size_t LinkedCellContainerDirect::findCellIndex(R3 vec) const {
    if (!fitsContainer(vec)) {
        SPDLOG_INFO("Position not not in container!");
        return cells.size();
    }
    if (!fitsDomain(vec)) {
        SPDLOG_INFO("Position of ghost particle!");
    }
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

std::vector<CellDirect*> LinkedCellContainerDirect::findAdjacentCells(size_t cell_idx) {
    std::vector<CellDirect*> adjacent_cells;
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
std::vector<const CellDirect*> LinkedCellContainerDirect::findAdjacentCells(size_t cell_idx) const {
    std::vector<const CellDirect*> adjacent_cells;
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

void LinkedCellContainerDirect::findBoundaryCells(const BoundaryType type,  // NOLINT
                                                  std::vector<const CellDirect*>& boundary_cells, size_t offset) const {
    switch (type) {
        case BoundaryType::UPPER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx =
                        ((num_cells[2] - (1U + offset)) * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::LOWER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (offset * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::LEFT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + offset;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::RIGHT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx =
                        (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + (num_cells[0] - (1U + offset));
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::FRONT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (offset * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::BACK: {
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
void LinkedCellContainerDirect::findBoundaryCells(const BoundaryType type,  // NOLINT
                                                  std::vector<CellDirect*>& boundary_cells, size_t offset) {
    switch (type) {
        case BoundaryType::UPPER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx =
                        ((num_cells[2] - (1U + offset)) * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::LOWER: {
            for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (offset * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::LEFT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + offset;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::RIGHT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t y = offset; y < num_cells[1] - offset; ++y) {
                    size_t idx =
                        (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + (num_cells[0] - (1U + offset));
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::FRONT: {
            for (size_t z = offset; z < num_cells[2] - offset; ++z) {
                for (size_t x = offset; x < num_cells[0] - offset; ++x) {
                    size_t idx = (z * num_cells[1] * num_cells[0]) + (offset * num_cells[0]) + x;
                    boundary_cells.push_back(&cells[idx]);
                }
            }
            return;
        }
        case BoundaryType::BACK: {
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

bool LinkedCellContainerDirect::fitsDomain(R3 v) const {
    return (v[0] >= 0.0 && v[0] <= domain_size[0]) && (v[1] >= 0.0 && v[1] <= domain_size[1]) &&
           (v[2] >= 0.0 && v[2] <= domain_size[2]);
}

bool LinkedCellContainerDirect::fitsContainer(R3 v) const {
    return (v[0] >= -cell_length[0] && v[0] <= domain_size[0] + cell_length[0]) &&
           (v[1] >= -cell_length[1] && v[1] <= domain_size[1] + cell_length[1]) &&
           (v[2] >= -cell_length[2] && v[2] <= domain_size[2] + cell_length[2]);
}

size_t LinkedCellContainerDirect::size() const {
    size_t n = 0;
    for (auto cell : cells) {
        n += cell.size();
    }
    return n;
}
bool LinkedCellContainerDirect::empty() const { return size() == 0; }

// modify
void LinkedCellContainerDirect::clear() {
    for (auto& cell : cells) {
        cell.clear();
    }
}

void LinkedCellContainerDirect::addParticle(Particle&& value) {
    if (!fitsContainer(value.getX())) {
        return;
    }
    cells[findCellIndex(value.getX())].addParticle(value);
}
void LinkedCellContainerDirect::addParticle(const Particle& value) {
    if (!fitsContainer(value.getX())) {
        return;
    }
    cells[findCellIndex(value.getX())].addParticle(value);
}

void LinkedCellContainerDirect::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    if (!fitsContainer(x_arg)) {
        return;
    }
    cells[findCellIndex(x_arg)].addParticle(Particle(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg));
}
void LinkedCellContainerDirect::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg,
                                            int type) {
    if (!fitsContainer(x_arg)) {
        return;
    }
    cells[findCellIndex(x_arg)].addParticle(Particle(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type));
}

void LinkedCellContainerDirect::eraseParticle(Particle* p) {
    size_t cell_idx = findCellIndex(p->getX());

    if (cell_idx < cells.size()) {
        auto it = std::find_if(cells[cell_idx].particles().begin(), cells[cell_idx].particles().end(),
                               [p](const Particle& particle) { return &particle == p; });  // NOLINT
        if (it != cells[cell_idx].particles().end()) {
            cells[cell_idx].removeParticle(it - cells[cell_idx].particles().begin());
        }
    }
}

void LinkedCellContainerDirect::updateParticlePosition(Particle& p, R3 new_x) {
    size_t old_cell_idx = findCellIndex(p.getX());
    if (!cells[old_cell_idx].fits(new_x)) {
        size_t new_cell_idx = findCellIndex(new_x);
        size_t old_idx = std::find(cells[old_cell_idx].particles().begin(), cells[old_cell_idx].particles().end(), p) -
                         cells[old_cell_idx].particles().begin();
        if (new_cell_idx == cells.size()) {
            cells[old_cell_idx].removeParticle(old_idx);
            return;
        }
        cells[new_cell_idx].addParticle(cells[old_cell_idx].removeParticle(old_idx));
    }

    p.getX() = new_x;
}

// normal iterators
LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::begin() {
    std::vector<CellDirect*> relevant_cells;
    relevant_cells.reserve(cells.size());
    for (auto& cell : cells) {
        relevant_cells.push_back(&cell);
    }
    return proximity_iterator{R3(), std::numeric_limits<double>::infinity(), cells.front().particles().begin(),
                              relevant_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::begin() const {
    std::vector<const CellDirect*> relevant_cells;
    relevant_cells.reserve(cells.size());
    for (const auto& cell : cells) {
        relevant_cells.push_back(&cell);
    }
    return const_proximity_iterator{R3(), std::numeric_limits<double>::infinity(), cells.front().particles().begin(),
                                    relevant_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::cbegin() const {
    std::vector<const CellDirect*> relevant_cells;
    relevant_cells.reserve(cells.size());
    for (const auto& cell : cells) {
        relevant_cells.push_back(&cell);
    }
    return const_proximity_iterator{R3(), std::numeric_limits<double>::infinity(), cells.front().particles().begin(),
                                    relevant_cells};
}
LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::end() {
    std::vector<CellDirect*> relevant_cells;
    relevant_cells.reserve(cells.size());
    for (auto& cell : cells) {
        relevant_cells.push_back(&cell);
    }
    return proximity_iterator{R3(), std::numeric_limits<double>::infinity(), cells.back().particles().end(),
                              relevant_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::end() const {
    std::vector<const CellDirect*> relevant_cells;
    relevant_cells.reserve(cells.size());
    for (const auto& cell : cells) {
        relevant_cells.push_back(&cell);
    }
    return const_proximity_iterator{R3(), std::numeric_limits<double>::infinity(), cells.back().particles().end(),
                                    relevant_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::cend() const {
    std::vector<const CellDirect*> relevant_cells;
    relevant_cells.reserve(cells.size());
    for (const auto& cell : cells) {
        relevant_cells.push_back(&cell);
    }
    return const_proximity_iterator{R3(), std::numeric_limits<double>::infinity(), cells.back().particles().end(),
                                    relevant_cells};
}

// proximity iterators
LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::proximityBegin(
    R3 center, double radius, [[maybe_unused]] size_t offset) {
    std::vector<CellDirect*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<CellDirect*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (CellDirect* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    return proximity_iterator{center, radius, nonempty_adjacent_cells.front()->particles().begin(),
                              nonempty_adjacent_cells};
}

LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::proximityEnd(R3 center, double radius) {
    std::vector<CellDirect*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<CellDirect*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (CellDirect* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    return proximity_iterator{center, radius, nonempty_adjacent_cells.back()->particles().end(),
                              nonempty_adjacent_cells};
}

LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::proximityBegin(
    R3 center, double radius, [[maybe_unused]] size_t offset) const {
    std::vector<const CellDirect*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<const CellDirect*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (const CellDirect* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    return const_proximity_iterator{center, radius, nonempty_adjacent_cells.front()->particles().begin(),
                                    nonempty_adjacent_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::proximityEnd(R3 center,
                                                                                            double radius) const {
    std::vector<const CellDirect*> adjacent_cells = findAdjacentCells(findCellIndex(center));
    std::vector<const CellDirect*> nonempty_adjacent_cells;

    nonempty_adjacent_cells.reserve(cells.size());
    for (const CellDirect* c : adjacent_cells) {
        if (!c->particles().empty()) {
            nonempty_adjacent_cells.push_back(c);
        }
    }

    return const_proximity_iterator{center, radius, nonempty_adjacent_cells.back()->particles().end(),
                                    nonempty_adjacent_cells};
}

// boundary and halo iterators
LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::haloBegin(
    const std::set<BoundaryType>& boundary_types) {
    std::vector<CellDirect*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::set<CellDirect*> seen;
    std::vector<CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.front()->particles().begin(), unique_and_nonempty_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::haloBegin(
    const std::set<BoundaryType>& boundary_types) const {
    std::vector<const CellDirect*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::set<const CellDirect*> seen;
    std::vector<const CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return const_proximity_iterator{
        R3{},
        std::numeric_limits<double>::infinity(),
        unique_and_nonempty_cells.front()->particles().begin(),
        unique_and_nonempty_cells,
    };
}
LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::haloEnd(
    const std::set<BoundaryType>& boundary_types) {
    std::vector<CellDirect*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::unordered_set<CellDirect*> seen;
    std::vector<CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.back()->particles().end(), unique_and_nonempty_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::haloEnd(
    const std::set<BoundaryType>& boundary_types) const {
    std::vector<const CellDirect*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 0);
    }

    std::unordered_set<const CellDirect*> seen;
    std::vector<const CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return const_proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                                    unique_and_nonempty_cells.back()->particles().end(), unique_and_nonempty_cells};
}

LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::boundaryBegin(
    const std::set<BoundaryType>& boundary_types) {
    std::vector<CellDirect*> relevant_cells;
    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<CellDirect*> seen;
    std::vector<CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.front()->particles().begin(), unique_and_nonempty_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::boundaryBegin(
    const std::set<BoundaryType>& boundary_types) const {
    std::vector<const CellDirect*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<const CellDirect*> seen;
    std::vector<const CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return const_proximity_iterator{
        R3{},
        std::numeric_limits<double>::infinity(),
        unique_and_nonempty_cells.front()->particles().begin(),
        unique_and_nonempty_cells,
    };
}
LinkedCellContainerDirect::proximity_iterator LinkedCellContainerDirect::boundaryEnd(
    const std::set<BoundaryType>& boundary_types) {
    std::vector<CellDirect*> relevant_cells;
    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<CellDirect*> seen;
    std::vector<CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return proximity_iterator{R3{}, std::numeric_limits<double>::infinity(),
                              unique_and_nonempty_cells.back()->particles().end(), unique_and_nonempty_cells};
}
LinkedCellContainerDirect::const_proximity_iterator LinkedCellContainerDirect::boundaryEnd(
    const std::set<BoundaryType>& boundary_types) const {
    std::vector<const CellDirect*> relevant_cells;

    for (const auto type : boundary_types) {
        findBoundaryCells(type, relevant_cells, 1);
    }

    std::unordered_set<const CellDirect*> seen;
    std::vector<const CellDirect*> unique_and_nonempty_cells;

    unique_and_nonempty_cells.reserve(cells.size());
    for (const CellDirect* c : relevant_cells) {
        if (seen.insert(c).second && !c->particles().empty()) {
            unique_and_nonempty_cells.push_back(c);
        }
    }

    return const_proximity_iterator{
        R3{},
        std::numeric_limits<double>::infinity(),
        unique_and_nonempty_cells.back()->particles().end(),
        unique_and_nonempty_cells,
    };
}

bool LinkedCellContainerDirect::isOnBoundary(Particle& p) {
    return (cells[findCellIndex(p.getX())].getType() == CellType::BOUNDARY);
}

R3 LinkedCellContainerDirect::getDomainSize() { return domain_size; }
