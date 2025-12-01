#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstddef>
#include <execution>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
#include "particles/boundaries/Boundary.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/SimpleContainer.h"
#include "particles/container/domain/Domain.h"
#include "physics/ForceSource.h"
#include "utils/Settings.h"

/**
 * @namespace mol_sim
 * @brief Namespace encompassing the project.
 *
 * Avoids name clashes in the case, that this project would be included in other projects.
 */
namespace mol_sim {

/**
 * @brief Class for running Simulation.
 * This Class implements a Builder Pattern. To run a full simulation, specify all needed parameters in the constructor
 * and then execute run().
 * @tparam containerType Type of container used for this simulation. Templated to work with Concept.
 * @tparam forceType Type of force source used for this simulation.
 */
template <ParticleContainer containerType, ForceSource forceType>
class Simulation {
   private:
    /**
     * @brief Simulation domain containing boundary conditions.
     */
    Domain domain;

    /**
     * @brief Container of Particles to simulate.
     */
    containerType& particles;

    /**
     * @brief Force source for calculating particle interactions.
     */
    forceType& force_source;

    /**
     * @brief Time step of simulation.
     */
    double delta_t;

    /**
     * @brief Start time of simulation.
     */
    double start_time;

    /**
     * @brief End time of simulation.
     */
    double end_time;

    /**
     * @brief Frequency of output file writing.
     */
    size_t frequency;

    /**
     * @brief Base name for output files.
     */
    std::string base_name;

    /**
     * @brief Cutoff radius for particles in proximity.
     */
    double cutoff_radius = std::numeric_limits<double>::infinity();

   public:
    /**
     * @brief Construct a new Simulation object and prepare for run() call.
     *
     * @param particles Container of particles to be simulated.
     * @param force_source Force source to be used in the simulation.
     * @param settings Simulation parameters.
     */
    Simulation(containerType& particles, forceType& force_source, SettingsParam& settings)
        : domain(std::move(settings.domain.value())),
          particles(particles),
          force_source(force_source),
          delta_t(settings.delta_t.value()),
          start_time(settings.start_time.value()),
          end_time(settings.end_time.value()),
          frequency(settings.frequency.value()),
          base_name(settings.base_name.value()),
          cutoff_radius(settings.cutoff.value()) {}

    /**
     * @brief Applies boundary conditions by creating ghost particles for reflecting boundaries.
     */
    void applyReflectingBoundaries() {
        std::vector<Particle> ghosts;
        // TODO: Optimze using boundary iterator
        auto it = particles.boundaryBegin();
        const auto end = particles.boundaryEnd();
        for (; it != end; ++it) {
            if (it->getType() != -1) {  // Don't create ghosts for ghost particles
                auto new_ghosts = domain.applyBoundary(*it);
                for (auto& ghost : new_ghosts) {
                    ghosts.push_back(std::move(ghost));
                }
            }
        }
        for (auto& ghost : ghosts) {
            particles.addParticle(std::move(ghost));
        }
    }

    /**
     * @brief Removes all particles in the Halo from the container.
     */
    void removeParticles() {
        // Collect indices of particles to remove using halo iterator
        std::vector<size_t> to_remove;
        for (auto it = particles.haloBegin(); it != particles.haloEnd(); ++it) {
            size_t idx = &(*it) - &particles[0];
            to_remove.push_back(idx);
        }

        // Sort in descending order to remove from end first (avoids index shifting issues)
        std::sort(to_remove.begin(), to_remove.end(), std::greater<size_t>());  // NOLINT

        // Remove particles using the standard vector iterator version
        for (size_t idx : to_remove) {
            particles.eraseParticle(particles.begin() + static_cast<std::ptrdiff_t>(idx));
        }
    }

    /**
     * @brief Calculates the forces of every particle for the next time step.
     */
    void calculateF() {
        for (auto& p : particles) {
            p.getOldF() = p.getF();
            p.getF() = Vector<double, 3>();
        }

        size_t idx = 1;
        for (auto it = particles.begin(); it != particles.end(); ++it, idx++) {
            Particle& p1 = *it;
            for (auto it_prox = particles.proximityBegin(p1.getX(), cutoff_radius, idx);
                 it_prox != particles.proximityEnd(p1.getX(), cutoff_radius); ++it_prox) {
                Particle& p2 = *it_prox;
                Vector<double, 3> force = force_source.applyForce(p1, p2);
                // Apply force directly (Newton's 3rd law: equal and opposite)
                p1.getF() = p1.getF() + force;
                p2.getF() = p2.getF() - force;
            }
        }
    }

    /**
     * @brief Calculates the positions of every particle for the next time step.
     */
    void calculateX() {
        for (auto it = particles.begin(); it != particles.end();) {
            const auto new_position =
                (*it).getX() + (delta_t * (*it).getV()) + ((0.5 * delta_t * delta_t / (*it).getM()) * (*it).getF());
            it = particles.updateParticlePosition(it, new_position);
        }
    }

    /**
     * @brief Calculates the velocities of every particle for the next time step.
     */
    void calculateV() {
        for (auto& p : particles) {
            p.getV() = p.getV() + ((0.5 * delta_t / p.getM()) * (p.getOldF() + p.getF()));
        }
    }

    /**
     * @brief Performs a full simulation run.
     */
    void run() {
        double current_time = start_time;
        [[maybe_unused]] int iteration = 0;

        // for this loop, we assume: current x, current f and current v are known
        while (current_time < end_time) {
            // 1. Calculate new positions
            calculateX();

            // 2. Remove OOB particles
            removeParticles();

            // 3. Apply reflecting boundaries (create ghost particles)
            applyReflectingBoundaries();

            // 4. Calculate forces (including ghost interactions)
            calculateF();

            // 5. Remove Halo particles
            removeParticles();

            // 6. Calculate new velocities
            calculateV();

            iteration++;
#ifndef DISABLE_IO
            if (iteration % frequency == 0) {
                try {
                    std::string out_name = base_name;
#ifdef ENABLE_VTK_OUTPUT
                    out_name += "_vtk";
                    VTKWriter writer;
#else
                    out_name += "_xyz";
                    XYZWriter writer;
#endif
                    writer.plotParticles(particles, out_name, iteration);
                } catch (...) {
                    SPDLOG_ERROR("Something went wrong with plotting the Particles.");
                }
            }
#endif
            SPDLOG_INFO("Iteration {} finished.", iteration);
            current_time += delta_t;
        }
        SPDLOG_INFO("Output written. Terminating...");
    }
};

}  // namespace mol_sim

#endif
