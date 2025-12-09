#include "particles/boundaries/Periodic.h"

namespace mol_sim {
std::array<short, 6> Periodic::corners;

Periodic::Periodic(BoundaryLocation location, R3 domain_size, R3 cell_size) noexcept
    : Boundary(location, BoundaryType::PERIODIC, domain_size), cell_size(cell_size) 
{
        std::array<short, 4> corners_to_mirror = {0};
        switch (location) {
            case BoundaryLocation::LEFT: 
                corners_to_mirror = {0, 2, 4, 6};
                break;
            case BoundaryLocation::RIGHT: 
                corners_to_mirror = {1, 3, 5, 7};
                break;
            case BoundaryLocation::UPPER: 
                corners_to_mirror = {4, 5, 6, 7};
                break;
            case BoundaryLocation::LOWER: 
                corners_to_mirror = {0, 1, 2, 3};
                break;
            case BoundaryLocation::FRONT: 
                corners_to_mirror = {0, 1, 4, 5};
                break;
            case BoundaryLocation::BACK: 
                corners_to_mirror = {2, 3, 6, 7};
                break;
            default:
                SPDLOG_ERROR("Unrecognized boundary location!");
                break;
        }
        for (size_t i = 0; i < 4; i++) {
            corners[corners_to_mirror[i]] = 1;
        }
}

std::optional<std::vector<Particle>> Periodic::applyBoundary(Particle& p, [[maybe_unused]] const ForceSource& force) noexcept {
    size_t axis = getAxis();
    int sign = getSign();
    std::vector<Particle> mirrored_particles;
   
    // 1) teleport OOB particles
    if (sign < 0 && p.getX()[axis] < 0) {
        p.getX()[axis] += domain_size[axis];
    } else if (sign > 0 && p.getX()[axis] > domain_size[axis]) {
        p.getX()[axis] -= domain_size[axis];
    }

    // 2) mirror the boundary particles (make sure to give them type 1)
    // 2.1) mirror the boundaries
    if (sign < 0 && p.getX()[axis] <= cell_size[axis] && p.getX()[axis] > 0) {
        Particle p_prime(p);
        p_prime.getX() = p.getX()[axis] + domain_size[axis];
        p_prime.getType() = 1;
        mirrored_particles.push_back(p_prime);   
    } else if (sign > 0 && p.getX()[axis] >= domain_size[axis] - cell_size[axis] 
                && p.getX()[axis] < domain_size[axis]) {
        Particle p_prime(p);
        p_prime.getX() = p.getX()[axis] - domain_size[axis];
        p_prime.getType() = 1;
        mirrored_particles.push_back(p_prime);   
    }

    // 2.2) mirror the corners
    for (size_t i = 0; i < 6; i++) {
        if (corners[i] == -1 * mark) {
            if (isInCorner(p)) {
                Particle p_prime(p);
                p_prime.getType() = 1;
                p_prime.getX() = p.getX() + getShift(i);
                mirrored_particles.push_back(p_prime);   
                corners[i] = mark;
            }
        }
    }
    mark *= -1;
    
    return mirrored_particles;
}

bool Periodic::isInCorner(Particle& p) const noexcept {
    R3 x = p.getX();
    return 
        ((x[0] >= 0 && x[0] <= cell_size[0]) || (x[0] <= domain_size[0] && x[0] >= domain_size[0] - cell_size[0])) && 
        ((x[1] >= 0 && x[1] <= cell_size[1]) || (x[1] <= domain_size[1] && x[1] >= domain_size[1] - cell_size[1])) && 
        ((x[2] >= 0 && x[2] <= cell_size[2]) || (x[2] <= domain_size[2] && x[2] >= domain_size[2] - cell_size[2])); 

}

R3 Periodic::getShift(size_t corner_idx) const noexcept {
    R3 shift = {0., 0., 0.};
    switch(corners[corner_idx]) {
        case 0:
            shift = domain_size;
            break; 
        case 1:
            shift = {-domain_size[0], domain_size[1], domain_size[2]};
            break; 
        case 2:
            shift = {domain_size[0], -domain_size[1], domain_size[2]};
            break; 
        case 3:
            shift = {-domain_size[0], -domain_size[1], domain_size[2]};
            break; 
        case 4:
            shift = {domain_size[0], domain_size[1], -domain_size[2]};
            break; 
        case 5:
            shift = {-domain_size[0], domain_size[1], -domain_size[2]};
            break; 
        case 6:
            shift = {domain_size[0], -domain_size[1], -domain_size[2]};
            break; 
        case 7:
            shift = -1 * domain_size;
            break;
        default:
            SPDLOG_ERROR("Unrecognized corner index!");
            break; 
    }
    return shift;
}

}  // namespace mol_sim
