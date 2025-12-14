#include "particles/boundaries/Periodic.h"

#include "particles/ParticleContainer.h"

namespace mol_sim {

Periodic::Periodic(BoundaryLocation location, R3 domain_size, double cutoff) noexcept
    : Boundary(location, BoundaryType::PERIODIC, domain_size), cutoff(cutoff) {
    /**
     * TODO: Optimization possible
     * 
     * Doing it this way introduces additional performance overhead
     * if used with a SimpleContainer if the domain_size is not divisible 
     * by the cutoff radius in at least one dimension.
     */
    std::array<size_t, 3> unused;
    LinkedCellContainer::computeCellsOrCorners(unused, halo_dimension, domain_size, cutoff);
}

bool Periodic::isInHalo(R3 x) const noexcept {
    LinkedCellContainer checker(domain_size, cutoff);
    return checker.fitsContainer(x) && !checker.fitsDomain(x);
}

bool Periodic::isOnBoundary(R3 x) const noexcept {
    size_t axis = getAxis();
    return ((x[axis] >= 0 && x[axis] <= halo_dimension[0]) ||
            (x[axis] <= domain_size[axis] && x[axis] >= domain_size[axis] - halo_dimension[axis])) &&
           (x[(axis + 1) % 3] >= 0 && x[(axis + 1) % 3] <= domain_size[(axis + 1) % 3]) &&
           (x[(axis + 2) % 3] >= 0 && x[(axis + 2) % 3] <= domain_size[(axis + 2) % 3]);
}

std::optional<std::vector<Particle>> Periodic::applyBoundary(  // NOLINT
    Particle& p, [[maybe_unused]] const ForceSource& force) noexcept {
    
    if (!isOnBoundary(p.getX())) { //only mirror particles in boundary (-> see task description)
        return std::nullopt;
    }

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
    /**
     * TODO: Optimization possible
     * e.g. via checking if domain_size > haloDimension (rules out half 
     * of the possible locations)
     */
    R3 offset = -1 * domain_size;
    for (size_t i = 0; i < 27; i++) {
        if (i == 13) {
            continue; //this is where p is right now.
        }
        if (i % 3 == 0 && i > 0) {
            offset[0] = offset[0] == 3 * domain_size[0] ? -domain_size[0] : offset[0] + domain_size[0];
            offset[2] = -domain_size[2];
        }
        if (i % 9 == 0 && i > 0) {
            offset[1] = -domain_size[1];
        }

        if (isInHalo(p.getX() + offset) && (p.getMirrorLocations() & 1 << i) == 0) {
            Particle p_prime(p);
            p_prime.getX() = p.getX() + offset;
            p_prime.getType() = 1;
            mirrored_particles.push_back(p_prime);
        }

        p.getMirrorLocations() |= 1 << i;
        offset[2] = offset[2] == 3 * domain_size[2] ? -domain_size[2] : offset[2] + domain_size[2];
    }  
    
    return mirrored_particles;
}

}  // namespace mol_sim
