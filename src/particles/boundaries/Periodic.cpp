#include "particles/boundaries/Periodic.h"

#include "particles/ParticleContainer.h"

namespace mol_sim {

Periodic::Periodic(BoundaryLocation location, R3 domain_size, double cutoff, size_t dimensions) noexcept
    : Boundary(location, BoundaryType::PERIODIC, domain_size), dimensions(dimensions) {
    /**
     * TODO: Optimization possible
     *
     * Doing it this way introduces additional performance overhead
     * if used with a SimpleContainer if the domain_size is not divisible
     * by the cutoff radius in at least one dimension.
     */
    std::array<size_t, 3> unused;
    LinkedCellContainer::computeCellsOrCorners(unused, halo_dimension, domain_size, cutoff);
    setShiftLookup(domain_size);
}

// NOLINTBEGIN
void Periodic::setShiftLookup(R3 domain_size) {
    // z inc -> x inc -> y inc
    // first x/z plane-------------------------------------------------------------
    shift_lookup[0] = {-domain_size[0], -domain_size[1], -domain_size[2]};  // 0
    shift_lookup[1] = {-domain_size[0], -domain_size[1], .0};               // 1
    shift_lookup[2] = {-domain_size[0], -domain_size[1], domain_size[2]};   // 2
    shift_lookup[3] = {.0, -domain_size[1], -domain_size[2]};               // 3
    shift_lookup[4] = {.0, -domain_size[1], .0};                            // 4
    shift_lookup[5] = {.0, -domain_size[1], domain_size[2]};                // 5
    shift_lookup[6] = {domain_size[0], -domain_size[1], -domain_size[2]};   // 6
    shift_lookup[7] = {domain_size[0], -domain_size[1], .0};                // 7
    shift_lookup[8] = {domain_size[0], -domain_size[1], domain_size[2]};    // 8
    // second x/z plane-----------------------------------------------------------
    shift_lookup[9] = {-domain_size[0], .0, -domain_size[2]};  // 9
    shift_lookup[10] = {-domain_size[0], .0, .0};              // 10
    shift_lookup[11] = {-domain_size[0], .0, domain_size[2]};  // 11
    shift_lookup[12] = {.0, .0, -domain_size[2]};              // 12
    shift_lookup[13] = {.0, .0, .0};                           // 13
    shift_lookup[14] = {.0, .0, domain_size[2]};               // 14
    shift_lookup[15] = {domain_size[0], .0, -domain_size[2]};  // 15
    shift_lookup[16] = {domain_size[0], .0, .0};               // 16
    shift_lookup[17] = {domain_size[0], .0, domain_size[2]};   // 17
    // third x/z plane------------------------------------------------------------
    shift_lookup[18] = {-domain_size[0], domain_size[1], -domain_size[2]};  // 18
    shift_lookup[19] = {-domain_size[0], domain_size[1], .0};               // 19
    shift_lookup[20] = {-domain_size[0], domain_size[1], domain_size[2]};   // 20
    shift_lookup[21] = {.0, domain_size[1], -domain_size[2]};               // 21
    shift_lookup[22] = {.0, domain_size[1], .0};                            // 22
    shift_lookup[23] = {.0, domain_size[1], domain_size[2]};                // 23
    shift_lookup[24] = {domain_size[0], domain_size[1], -domain_size[2]};   // 24
    shift_lookup[25] = {domain_size[0], domain_size[1], .0};                // 25
    shift_lookup[26] = {domain_size[0], domain_size[1], domain_size[2]};    // 26
}
// NOLINTEND

/**
 * TODO: Edge case:
 * Imagine a situation where we have say a 2D simulation where LEFT and RIGHT are periodic
 * and the other boundaries aren't. Inside the iteration in Domain.cpp:applyBoundary either the left
 * boundary or the right boundary will execute their applyBoundary() before the other. I.e. either
 * its left->applyBoundary then right->applyBoundary or vice versa. Let's say the order is LEFT->RIGHT.
 * Particles that are teleported from the left halo to the right boundary will be mirrored (if necessary)
 * by the right boundary in the subsequent call to right->applyBoundary in Domain.cpp:applyBoundary.
 * However, the same does NOT hold for the right boundary. That is: in general particles that are teleported
 * from the right halo to the left boundary will NOT be mirrored in that same timestep (unless the particle
 * is in an edge or a corner). They will be mirrored (if necessary) in the next timestep.
 * This *could* lead to instabilitiees in the simulation. Assuming reasonably sized simulation parameters
 * in particular delta_t, sigma and epsilon of the relevant particles, the missing of these mirror particles
 * in the one timestep may not be critical. In the interest of not entangling this code even more
 * (we would need to pass in a SettingsParam& to the constructor here + have to check for the
 * types of the other boundaries) and not overcomplicating this code (unless really necessary)
 * I will leave this behavior in for now.
 */
std::optional<std::vector<Particle>> Periodic::applyBoundary(  // NOLINT
    Particle& p, [[maybe_unused]] const PairwiseForceSource& force) noexcept {
    teleportParticleIfOOB(p);
    if (!isOnBoundary(p.getX(), getAxis(), getSign())) {  // only mirror particles in correct boundary region
        return std::nullopt;
    }
    auto mirrored_particles = mirrorParticle(p);
    if (mirrored_particles.empty()) {
        return std::nullopt;
    }
    return mirrored_particles;
}

void Periodic::teleportParticleIfOOB(Particle& p) {
    R3 position_shift = {.0, .0, .0};
    size_t axis = getAxis();
    int sign = getSign();
    for (size_t i = 0; i < dimensions; i++, axis = (axis + 1) % 3) {
        if ((i == 0 ? sign < 0 : true) && p.getX()[axis] < 0) {
            position_shift[axis] = domain_size[axis];
        } else if ((i == 0 ? sign > 0 : true) && p.getX()[axis] > domain_size[axis]) {
            position_shift[axis] = -domain_size[axis];
        }
        if (i == 0 && position_shift[axis] == 0) {
            break;  // if the particle is not in halo region of the periodic boundary, definitely DONT shift!
        }
    }
    p.getX() = p.getX() + position_shift;
}

/**
 * TODO: make this work with an std::array<Particle, 4>
 */
std::vector<Particle> Periodic::mirrorParticle(Particle& p) {
    std::vector<Particle> mirrored_particles;
    size_t location_idx = getBoundaryLocationIdx();

    // 1) definitely mirror particle to other side once
    R3 shift = {.0, .0, .0};
    size_t axis = getAxis();
    shift[axis] = -getSign() * domain_size[axis];
    addMirrorParticle(p.getX() + shift, MIRROR_IDX_LOOKUP[(dimensions == 2 ? 58 : 48) + location_idx], p,
                      mirrored_particles);

    // 2 handle edges and corners
    size_t idx = getIdx(p);
    if (idx <= 3) {  // if p in some 3D edge = 2D corner
        uint8_t mirror_idx = MIRROR_IDX_LOOKUP[idx + (dimensions == 2 ? 54 : 8 * location_idx)];
        addMirrorParticle(p.getX() + shift_lookup[mirror_idx], mirror_idx, p, mirrored_particles);
    } else if (dimensions == 3 && 4 <= idx && idx <= 7) {  // if p in some 3D corner
        uint8_t mirror_idx1 = MIRROR_IDX_LOOKUP[(8 * location_idx) + (idx % 4)];
        uint8_t mirror_idx2 = MIRROR_IDX_LOOKUP[(8 * location_idx) + ((idx + 1) % 4)];
        uint8_t mirror_idx = MIRROR_IDX_LOOKUP[idx + (8 * location_idx)];
        addMirrorParticle(p.getX() + shift_lookup[mirror_idx1], mirror_idx1, p, mirrored_particles);
        addMirrorParticle(p.getX() + shift_lookup[mirror_idx2], mirror_idx2, p, mirrored_particles);
        addMirrorParticle(p.getX() + shift_lookup[mirror_idx], mirror_idx, p, mirrored_particles);
    }
    return mirrored_particles;
}

void Periodic::addMirrorParticle(const R3& mirrorLocation, size_t mirrorIdx, Particle& p,
                                 std::vector<Particle>& mirrored_particles) {
    if ((p.getMirrorLocations() & (1 << mirrorIdx)) == 0) {
        Particle mirror_particle(p);
        mirror_particle.getX() = mirrorLocation;
        mirror_particle.getType() = 1;
        p.getMirrorLocations() |= (1 << mirrorIdx);
        mirrored_particles.push_back(mirror_particle);
    }
}

// NOLINTBEGIN
size_t Periodic::getIdx(Particle& p) {
    size_t axis = getAxis();
    if (dimensions == 3) {  // handle 3D corner and edge indicies
        size_t axis1 = 0;
        size_t axis2 = 0;
        switch (axis) {
            case 0: {
                axis1 = 1;
                axis2 = 2;
                break;
            }
            case 1: {
                axis1 = 0;
                axis2 = 2;
                break;
            }
            case 2: {
                axis1 = 0;
                axis2 = 1;
                break;
            }
            default: {
                SPDLOG_ERROR("Unknown axis!");
            }
        }
        bool axis1_small = p.getX()[axis1] <= halo_dimension[axis1];
        bool axis1_big = p.getX()[axis1] >= domain_size[axis1] - halo_dimension[axis1];
        bool axis2_small = p.getX()[axis2] <= halo_dimension[axis2];
        bool axis2_big = p.getX()[axis2] >= domain_size[axis2] - halo_dimension[axis2];
        if (axis1_small && axis2_small) return 4;
        if (axis1_big && axis2_small) return 5;
        if (axis1_big && axis2_big) return 6;
        if (axis1_small && axis2_big) return 7;
        if (axis1_small) return 0;
        if (axis2_small) return 1;
        if (axis1_big) return 2;
        if (axis2_big) return 3;
    } else if (dimensions == 2) {  // handle 2D corner indicies
        bool x_small = p.getX()[0] <= halo_dimension[0] && p.getX()[0] >= 0;
        bool x_big = p.getX()[0] >= domain_size[0] - halo_dimension[0] && p.getX()[0] <= domain_size[0];
        bool y_small = p.getX()[1] <= halo_dimension[1] && p.getX()[1] >= 0;
        bool y_big = p.getX()[1] >= domain_size[1] - halo_dimension[1] && p.getX()[1] <= domain_size[1];
        if (x_small && y_small) return 0;
        if (x_big && y_small) return 1;
        if (x_big && y_big) return 2;
        if (x_small && y_big) return 3;
    }
    return 8;
}
// NOLINTEND

bool Periodic::isOnBoundary(R3 x, size_t axis, int sign) const noexcept {
    bool result = true;
    for (size_t i = 0; i < dimensions; i++, axis = (axis + 1) % 3) {
        if (i == 0) {
            result &= (sign < 0)
                          ? (x[axis] >= 0 && x[axis] <= halo_dimension[axis])
                          : (x[axis] <= domain_size[axis] && x[axis] >= domain_size[axis] - halo_dimension[axis]);
        } else {
            result &= (x[axis] >= 0 && x[axis] <= domain_size[axis]);
        }
    }
    return result;
}
}  // namespace mol_sim
