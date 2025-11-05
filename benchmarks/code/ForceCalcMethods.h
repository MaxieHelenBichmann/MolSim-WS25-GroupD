#ifndef FORCE_CALC_H
#define FORCE_CALC_H
#include <memory>
#include <vector>

#include "particles/ParticleContainer.h"
#include "physics/ForceSource.h"
namespace mol_sim {
template <ParticleContainer containerType>
void calculateF(containerType& particles, std::unique_ptr<ForceSource>& force_source) {
    std::vector<Vector<double, 3>> forces(particles.size(), Vector<double, 3>());

    for (auto& p : particles) {
        p.getOldF() = p.getF();
        p.getF() = Vector<double, 3>();
    }

    for (size_t i = 0; i < particles.size(); i++) {
        Particle& p1 = particles[i];
        // compute row
        for (size_t j = i + 1; j < particles.size(); j++) {
            Particle& p2 = particles[j];
            Vector<double, 3> force = force_source->calculateForce(p1, p2);
            forces[j] = force;
            p1.getF() = p1.getF() + force;
        }
        // apply row in column
        for (size_t j = i + 1; j < particles.size(); j++) {
            particles[j].getF() = particles[j].getF() - forces[j];
        }
    }
}

/**
 * @brief Calculates the forces of every particle for the next time step.
 * Calculates the forces of every particle. for the next time step. Using the specified force source and delta_t.
 *
 * ONLY HERE FOR BENCHMARKING REASONS.
 */
template <ParticleContainer containerType>
void calculateFUnoptimized(containerType& particles, std::unique_ptr<ForceSource>& force_source) {
    for (auto& p1 : particles) {
        p1.getOldF() = p1.getF();
        p1.getF() = Vector<double, 3>();
        for (auto& p2 : particles) {
            if (p1 == p2) {
                continue;
            }
            p1.getF() = p1.getF() + force_source->calculateForce(p1, p2);
        }
    }
}
}  // namespace mol_sim
#endif
