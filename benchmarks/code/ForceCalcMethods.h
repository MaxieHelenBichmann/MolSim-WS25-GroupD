#ifndef FORCE_CALC_H
#define FORCE_CALC_H
#include <memory>
#include <vector>

#include "AbstractForce.h"
#include "particles/ParticleContainer.h"
#include "physics/ForceSource.h"
namespace mol_sim {

/**
 * @deprecated ONLY USED FOR BENCHMARKING
 * @brief Copy of the calculateF Method in Simulate.
 *
 * @tparam containerType Type of container used.
 * @param particles Particle Container to calculate forces for.
 * @param force_source Force source to be used for calculation
 */
template <ParticleContainer containerType>
void calculateF(containerType& particles, std::unique_ptr<ForceAbstract>& force_source) {
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
            Vector<double, 3> force = force_source->applyForce(p1, p2);
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
 * @deprecated ONLY USED FOR BENCHMARKING
 * @brief Old unoptimized version of the calculateF Method in Simulation.
 *
 * @tparam containerType Type of container used.
 * @param particles Particle Container to calculate forces for.
 * @param force_source Force source to be used for calculation
 */
template <ParticleContainer containerType>
void calculateFUnoptimized(containerType& particles, std::unique_ptr<ForceAbstract>& force_source) {
    for (auto& p1 : particles) {
        p1.getOldF() = p1.getF();
        p1.getF() = Vector<double, 3>();
        for (auto& p2 : particles) {
            if (p1 == p2) {
                continue;
            }
            p1.getF() = p1.getF() + force_source->applyForce(p1, p2);
        }
    }
}
}  // namespace mol_sim
#endif
