#include "PeriodicOld.h"

namespace mol_sim {

std::unique_ptr<LinkedCellContainer> PeriodicOld::checker;

PeriodicOld::PeriodicOld(BoundaryLocation location, R3 domain_size, double cutoff, size_t dimensions) noexcept
    : Boundary(location, BoundaryType::PERIODIC, domain_size) {
    /**
     * TODO: Optimization possible
     *
     * Doing it this way introduces additional performance overhead
     * if used with a SimpleContainer if the domain_size is not divisible
     * by the cutoff radius in at least one dimension.
     */
    if (checker == nullptr) {
        checker = std::make_unique<LinkedCellContainer>(domain_size, cutoff);
    }
    std::array<size_t, 3> unused;
    LinkedCellContainer::computeCellsOrCorners(unused, halo_dimension, domain_size, cutoff);
    is2D = dimensions == 2;
}

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
 * in the one timestep may not be critical. In th interest of not entangling this code even more
 * (we would need to pass in a SettingsParam& to the constructor here + have to check for the
 * types of the other boundaries) and not overcomplicating this code (unless really necessary)
 * I will leave this behavior in for now.
 */
void PeriodicOld::applyBoundary(  // NOLINT
    Particle& p, [[maybe_unused]] const ForceSource& force) noexcept {
    teleportParticleIfOOB(p);
    if (!isOnBoundary(p.getX(), getAxis(), getSign())) {  // only mirror particles in correct boundary region
        return;
    }
    mirrorParticle(p);
}

void PeriodicOld::teleportParticleIfOOB(Particle& p) {
    R3 position_shift = {.0, .0, .0};
    size_t axis = getAxis();
    int sign = getSign();
    for (size_t i = 0; is2D ? i < 2 : i < 3; i++, axis = (axis + 1) % 3) {
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
 * TODO: Optimization possible
 * e.g. via checking if domain_size > haloDimension (rules out half
 * of the possible locations)
 * Could also reserve vector or not use a vector for memory efficiency
 */
void PeriodicOld::mirrorParticle(Particle& p) {
    size_t stride = is2D ? 3 : 1;
    size_t i = is2D ? 1 : 0;
    R3 offset = -1 * domain_size;
    if (is2D) {
        offset[2] = 0;
    }
    for (; i < 27; i += stride) {
        // add mirror particle if necessary
        if (i != 13 && isInHalo(p.getX() + offset) && (p.getMirrorLocations() & (1 << i)) == 0) {
            Particle p_prime(p);
            p.getMirrorPositions().emplace_back(p.getX() + offset);
            p.getMirrorLocations() |= (1 << i);
        }

        // update offset
        updateOffset(offset, i);
    }
}

bool PeriodicOld::isOnBoundary(R3 x, size_t axis, int sign) const noexcept {
    bool result = true;
    for (size_t i = 0; is2D ? i < 2 : i < 3; i++, axis = (axis + 1) % 3) {
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

bool PeriodicOld::isInHalo(R3 x) const noexcept {
    size_t axis = 0;
    int sign = 0;
    for (size_t i = 0; is2D ? i < 2 : i < 3; i++) {
        if (x[i] == 0) {
            axis = i;
            sign = -1;
            break;
        }
        if (x[i] == domain_size[i]) {
            axis = i;
            sign = 1;
            break;
        }
    }
    return (checker->fitsContainer(x) && (!checker->fitsDomain(x)))  // halo around container OR
           || (sign != 0 && isOnBoundary(x, axis, sign));            // the boundary itself
}
/**
 * TODO: Does this need to be recalculated or could we store a static array of all possible offsets to optimize this a
 * bit?
 */
void PeriodicOld::updateOffset(R3& offset, size_t i) {
    if (is2D) {
        if (offset[0] == domain_size[0]) {
            offset[0] = -domain_size[0];
            offset[1] += domain_size[1];
        } else {
            offset[0] += domain_size[0];
        }
    } else {
        if ((i + 1) % 3 == 0) {
            offset[0] = offset[0] == domain_size[0] ? -domain_size[0] : offset[0] + domain_size[0];
            offset[2] = -domain_size[2];
        }
        if ((i + 1) % 9 == 0) {
            offset[1] += domain_size[1];
        }
        if ((i + 1) % 3 != 0) {
            offset[2] += domain_size[2];
        }
    }
}

}  // namespace mol_sim
