#include "particles/container/LinkedCellContainer.h"

using namespace mol_sim;

// ------------------- Cell methods -------------------

LinkedCellContainer::Cell::Cell(CellType cell_type, std::array<double, 6> bounds) : type(cell_type), bounds(bounds) {}

void LinkedCellContainer::Cell::addParticle(Particle* value) { data.push_back(value); }
void LinkedCellContainer::Cell::removeParticle(size_t idx) { data.erase(data.begin() + idx); }  // NOLINT
void LinkedCellContainer::Cell::removeParticle(Particle* p) {
    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i] == p) {
            data.erase(data.begin() + i);  // NOLINT
            return;
        }
    }
}
void LinkedCellContainer::Cell::clear() { data.clear(); }
Particle* LinkedCellContainer::Cell::operator[](size_t idx) { return data[idx]; }
bool LinkedCellContainer::Cell::fits(Particle* p) const {
    return bounds[0] <= p->getX()[0] && p->getX()[0] < bounds[1] && bounds[2] <= p->getX()[1] &&
           p->getX()[1] < bounds[3] && bounds[4] <= p->getX()[2] && p->getX()[2] < bounds[5];
}
size_t LinkedCellContainer::Cell::size() { return data.size(); }
LinkedCellContainer::Cell::CellType LinkedCellContainer::Cell::getType() { return type; }

// ------------------- LinkedCellContainer methods -------------------

LinkedCellContainer::LinkedCellContainer(R3 domain_size, double cutoff_radius)
    : domain_size(domain_size), cutoff_radius(cutoff_radius) {
    num_cells[0] = static_cast<size_t>(std::ceil(domain_size[0] / cutoff_radius));
    num_cells[1] = static_cast<size_t>(std::ceil(domain_size[1] / cutoff_radius));
    num_cells[2] = static_cast<size_t>(std::ceil(domain_size[2] / cutoff_radius));

    cells.reserve(num_cells[0] * num_cells[1] * num_cells[2]);
    for (size_t z = 0; z < num_cells[2]; ++z) {
        for (size_t y = 0; y < num_cells[1]; ++y) {
            for (size_t x = 0; x < num_cells[0]; ++x) {
                std::array<double, 6> bounds = {static_cast<double>(x) * cutoff_radius,
                                                static_cast<double>((x + 1)) * cutoff_radius < domain_size[0]
                                                    ? static_cast<double>((x + 1)) * cutoff_radius
                                                    : domain_size[0],
                                                static_cast<double>(y) * cutoff_radius,
                                                static_cast<double>((y + 1)) * cutoff_radius < domain_size[1]
                                                    ? static_cast<double>((y + 1)) * cutoff_radius
                                                    : domain_size[1],
                                                static_cast<double>(z) * cutoff_radius,
                                                static_cast<double>((z + 1)) * cutoff_radius < domain_size[2]
                                                    ? static_cast<double>((z + 1)) * cutoff_radius
                                                    : domain_size[2]};
                cells.emplace_back(Cell::CellType::INNER, bounds);
            }
        }
    }
}

void LinkedCellContainer::switchCell(Particle& p, size_t old_cell_idx, size_t new_cell_idx) {
    cells[old_cell_idx].removeParticle(&p);
    cells[new_cell_idx].addParticle(&p);
}

void LinkedCellContainer::switchCell(Particle* p, size_t old_cell_idx, size_t new_cell_idx) {
    cells[old_cell_idx].removeParticle(p);
    cells[new_cell_idx].addParticle(p);
}

void LinkedCellContainer::switchCell(size_t p, size_t old_cell_idx, size_t new_cell_idx) {
    Particle* particle_ptr = cells[old_cell_idx][p];
    cells[old_cell_idx].removeParticle(particle_ptr);
    cells[new_cell_idx].addParticle(particle_ptr);
}

size_t LinkedCellContainer::findCellIndex(Particle* p) {
    size_t z = std::floor(p->getX()[2] / cutoff_radius);
    size_t y = std::floor(p->getX()[1] / cutoff_radius);
    size_t x = std::floor(p->getX()[0] / cutoff_radius);

    return (z * num_cells[1] * num_cells[0]) + (y * num_cells[0]) + x;
}

bool LinkedCellContainer::fits(R3 v) const {
    return v[0] <= domain_size[0] && v[1] <= domain_size[1] && v[2] <= domain_size[2];
}

Particle& LinkedCellContainer::operator[](size_t idx) { return data[idx]; }
const Particle& LinkedCellContainer::operator[](size_t idx) const { return data[idx]; }

size_t LinkedCellContainer::size() const { return data.size(); }
bool LinkedCellContainer::empty() const { return data.empty(); }

// modify
void LinkedCellContainer::clear() {
    data.clear();
    for (auto& cell : cells) {
        cell.clear();
    }
}
void LinkedCellContainer::reserve(size_t n) { data.reserve(n); }

void LinkedCellContainer::addParticle(Particle&& value) {
    if (!fits(value.getX())) {
        return;
    }
    data.push_back(value);
    Particle* p = &data.back();  // NOLINT
    cells[findCellIndex(p)].addParticle(p);
}
void LinkedCellContainer::addParticle(const Particle& value) {
    if (!fits(value.getX())) {
        return;
    }
    data.push_back(value);
    Particle* p = &data.back();  // NOLINT
    cells[findCellIndex(p)].addParticle(p);
}

void LinkedCellContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    if (!fits(x_arg)) {
        return;
    }
    data.emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
    Particle* p = &data.back();  // NOLINT
    cells[findCellIndex(p)].addParticle(p);
}
void LinkedCellContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg,
                                      int type) {
    if (!fits(x_arg) && type != -1) {
        return;
    }
    data.emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
    Particle* p = &data.back();  // NOLINT
    cells[findCellIndex(p)].addParticle(p);
}

// normal iterators
std::vector<Particle>::iterator LinkedCellContainer::begin() { return data.begin(); }
std::vector<Particle>::const_iterator LinkedCellContainer::begin() const { return data.begin(); }
std::vector<Particle>::const_iterator LinkedCellContainer::cbegin() const { return data.cbegin(); }
std::vector<Particle>::iterator LinkedCellContainer::end() { return data.end(); }
std::vector<Particle>::const_iterator LinkedCellContainer::end() const { return data.end(); }
std::vector<Particle>::const_iterator LinkedCellContainer::cend() const { return data.cend(); }

// proximity iterators
LinkedCellContainer::proximity_iterator LinkedCellContainer::proximityBegin(R3 center, double radius, size_t offset) {
    (void)center;
    (void)radius;
    (void)offset;
    return {};
}

LinkedCellContainer::proximity_iterator LinkedCellContainer::proximityEnd(R3 center, double radius) {
    (void)center;
    (void)radius;
    return {};
}

LinkedCellContainer::const_proximity_iterator LinkedCellContainer::proximityBegin(R3 center, double radius,
                                                                                  size_t offset) const {
    (void)center;
    (void)radius;
    (void)offset;
    return {};
}
LinkedCellContainer::const_proximity_iterator LinkedCellContainer::proximityEnd(R3 center, double radius) const {
    (void)center;
    (void)radius;
    return {};
}

// static_assert(ParticleContainer<LinkedCellContainer>);

bool LinkedCellContainer::isOnBoundary(Particle& p) {
    return (cells[findCellIndex(&p)].getType() == Cell::CellType::BOUNDARY);
}

void LinkedCellContainer::removeParticle(Particle& p) { 
    cells[findCellIndex(&p)].removeParticle(&p); 
}

R3 LinkedCellContainer::getDomainSize() {
    return domain_size;
}