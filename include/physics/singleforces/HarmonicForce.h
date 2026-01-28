#ifndef HARMONIC_FORCE_H
#define HARMONIC_FORCE_H

#include <cstddef>
#include <numbers>

#include "particles/Particle.h"
#include "particles/container/ContainerRef.h"
#include "physics/singleforces/SingleForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate harmonic spring forces for membrane particles.
 *
 * A SingleForceSource that calculates forces based on harmonic potential between neighboring particles in a membrane
 * structure. Only applies to membrane particles (types 2 and 4).
 */
class HarmonicForce : public SingleForceSource {
   private:
    /**
     * @brief Stiffness constant of the harmonic spring.
     */
    const double K;
    /**
     * @brief Equilibrium bond length between direct neighbors.
     */
    const double R_0;

    ContainerRef particles;

   public:
    /**
     * @brief Construct a new Harmonic Force calculator.
     *
     * @param k Stiffness constant.
     * @param r_0 Equilibrium bond length.
     */
    HarmonicForce(double k, double r_0) : K(k), R_0(r_0) {}
    /**
     * @brief Calculates the harmonic spring force acting on a membrane particle.
     *
     * @param p1 Particle whose force is to be calculated.
     * @return Vector<double, 3> Total harmonic force from all neighbors.
     *
     * Calculates forces from direct neighbors using distance R_0 and diagonal neighbors using sqrt(2)*R_0.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1) const noexcept override {
        Vector<double, 3> force = {0., 0., 0.};
        if (p1.getType() != 2 && p1.getType() != 4) {
            return force;
        }
        const auto& neighbors = p1.getNeighbors();
        // compute direct neighbor influence
        for (size_t i = 0; i < 4; i++) {
            const auto neighbor = neighbors[i];
            if (neighbor.has_value()) {
                const Particle& p2 = particles[neighbor.value()];
                const R3 diff = p2.getX() - p1.getX();
                const double dist = diff.euclidNorm();
                const double scalar = (K * 0.5 * (dist - R_0)) / dist;
                force += scalar * diff;
            }
        }
        // compute diagonal neighbor influence
        for (size_t i = 4; i < 8; i++) {
            const auto neighbor = neighbors[i];
            if (neighbor.has_value()) {
                const Particle& p2 = particles[neighbor.value()];
                const R3 diff = p2.getX() - p1.getX();
                const double dist = diff.euclidNorm();
                const double scalar = (K * 0.5 * (dist - std::numbers::sqrt2 * R_0)) / dist;
                force += scalar * diff;
            }
        }
        return force;
    }

    /**
     * @brief Sets the particle container reference for force calculations.
     *
     * @param container Reference to the particle container.
     */
    void setContainer(ContainerRef container) { particles = container; }

    /**
     * @brief Gets the type of this force source.
     *
     * @return SingleForce type identifier (HARMONIC).
     */
    [[nodiscard]] SingleForce getType() const override { return HARMONIC; }
};
}  // namespace mol_sim

#endif
