#ifndef VELOCITYREFLECT_H
#define VELOCITYREFLECT_H

#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "utils/Vector.h"

namespace mol_sim {
class VelocityReflect : public Boundary {
   public:
    VelocityReflect(BoundaryLocation location, R3 domain_size) noexcept;
    ~VelocityReflect() override = default;

    /**
     * @brief The applyBoundary routine of VelocityReflect. If a particle moved beyond the domain in the current
     * time-step, this boundary takes the vector starting at where the particle intersected the domain and ending
     * at the new position of the particle and flips it around so the particle stays within the domain.
     * In other words this boundary condition can be thought of as throwing a ball against a flat wall.
     *
     * @param p The particle we want to apply the boundary condition on.
     * @param force The force source used  in this boundary condition (irrelevant in this case).
     * @return std::nullopt always, since no new particles are generated here.
     */
    void applyBoundary(Particle& p, [[maybe_unused]] const PairwiseForceSource& force) noexcept override;
};
}  // namespace mol_sim
#endif
