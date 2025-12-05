#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "exceptions/SimulationException.h"
#include "io/OutputWriter.h"
#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
#include "particles/container/domain/Domain.h"
#include "physics/ForceSource.h"
#include "utils/Settings.h"
#include "utils/Vector.h"

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
template <ParticleContainer containerType>
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
    const ForceSource& force_source;
    /**
     * @brief Writer used for output.
     */
    const OutputWriter& writer;
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

    size_t dimensions;
    /**
     * @brief Cutoff radius for particles in proximity.
     */
    double cutoff_radius;

    double total_energy;

    double target_temp;
    double delta_temp;

   public:
    /**
     * @brief Construct a new Simulation object and prepare for run() call.
     *
     * @param particles Container of particles to be simulated.
     * @param force_source Force source to be used in the simulation.
     * @param settings Simulation parameters.
     */
    Simulation(containerType& particles, const ForceSource& force_source, SettingsParam& settings,
               const OutputWriter& writer)
        : domain(std::move(settings.domain)),
          particles(particles),
          force_source(force_source),
          writer(writer),
          delta_t(settings.delta_t),
          start_time(settings.start_time),
          end_time(settings.end_time),
          frequency(settings.frequency),
          base_name(settings.base_name),
          dimensions(settings.dimensions),
          cutoff_radius(settings.cutoff),
          target_temp(settings.target_temp),
          delta_temp(settings.delta_temp) {
        for (const Particle& p : particles) {
            total_energy += p.getM() * R3::scalarProduct(p.getV(), p.getV());
        }
        total_energy *= 0.5;
    }

    /**
     * @brief Removes all particles in the Halo from the container.
     */
    void removeParticles() {
        // Collect indices of particles to remove using halo iterator
        SPDLOG_DEBUG("Container has currently {} particles before erase", particles.size());
        std::vector<size_t> to_remove;
        for (auto it = particles.haloBegin(); it != particles.haloEnd(); ++it) {
            size_t idx = &(*it) - &particles[0];
            to_remove.push_back(idx);
        }

        SPDLOG_DEBUG("Added {} (ghost) particles to remove", to_remove.size());

        // Sort in descending order to remove from end first (avoids index shifting issues)
        std::sort(to_remove.begin(), to_remove.end(), std::greater<size_t>());  // NOLINT

        // Remove particles using the standard vector iterator version
        for (size_t idx : to_remove) {
            particles.eraseParticle(particles.begin() + static_cast<std::ptrdiff_t>(idx));
        }
        SPDLOG_DEBUG("Container has currently {} particles after erase", particles.size());
    }

    /**
     * @brief Applies the necessary boundary conditions to the particles.
     */
    void applyBoundaries() {
        for (auto it = particles.begin(); it != particles.end();) {
            (*it).getOldF() = (*it).getF();
            (*it).getF() = Vector<double, 3>();
            // TODO Optimization to only call this for relevant particles
            domain.applyBoundary(*it, force_source);
            // TODO: bit of an ugly workaround for now.
            R3 new_position = (*it).getX();
            (*it).getX() = (*it).getOldX();
            it = particles.updateParticlePosition(it, new_position);
        }
    }

    /**
     * @brief Calculates the forces of every particle for the next time step.
     */
    void calculateF() {
        size_t idx = 0;
        for (auto it = particles.begin(); it != particles.end(); ++it, idx++) {
            Particle& p1 = *it;
            for (auto it_prox = particles.proximityBegin(p1.getX(), idx); it_prox != particles.proximityEnd(p1.getX());
                 ++it_prox) {
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
        for (auto& p : particles) {
            p.getOldX() = p.getX();
            p.getX() = p.getX() + (delta_t * p.getV()) + ((0.5 * delta_t * delta_t / p.getM()) * p.getF());
        }
    }

    /**
     * @brief Calculates the velocities of every particle for the next time step.
     */
    void calculateV(double scalar_factor) {
        double curr_energy = 0;
        for (auto& p : particles) {
            R3 new_v = scalar_factor * (p.getV() + ((0.5 * delta_t / p.getM()) * (p.getOldF() + p.getF())));
            p.getV() = new_v;
            curr_energy += p.getM() * R3::scalarProduct(new_v, new_v);
        }
        total_energy = 0.5 * curr_energy;
    }

    /**
     * @brief Performs a full simulation run.
     * @throws SimulationException if an error occurs during output writing.
     */
    void run() {
        double current_time = start_time;
        [[maybe_unused]] int iteration = 0;

        SPDLOG_INFO("Starting simulation: {} particles, t=[{}, {}], dt={}", particles.size(), start_time, end_time,
                    delta_t);

        // for this loop, we assume: current x, current f and current v are known
        while (current_time < end_time) {
            SPDLOG_DEBUG("Iteration {}: Updating {} particle positions", iteration + 1, particles.size());
            // 1. Calculate new positions
            calculateX();

            // 2. Apply boundaries
            applyBoundaries();
            // 3. Remove OOB particles
            removeParticles();

            // 4. Calculate forces (including ghost interactions)
            SPDLOG_DEBUG("Iteration {}: Calculating forces for {} particles", iteration + 1, particles.size());
            calculateF();

            // 5. Calculate thermostat factor
            double curr_temp = (2.0 * total_energy) / (dimensions * particles.size());

            double clamped_target = curr_temp + std::clamp((target_temp - curr_temp), -delta_temp, delta_temp);
            double thermo_factor = sqrt(clamped_target / curr_temp);
            // 6. Calculate new velocities
            calculateV(thermo_factor);

            iteration++;
#ifndef DISABLE_IO
            if (iteration % frequency == 0) {
                try {
                    std::string out_name = base_name;
                    writer.plotParticles(particles, out_name, iteration);
                } catch (const std::runtime_error& e) {
                    SPDLOG_ERROR("Failed to plot particles at iteration {}: {}", iteration, e.what());
                    throw SimulationException("Error while plotting Particles: " + std::string(e.what()));
                }
            }
#endif
            SPDLOG_DEBUG("Iteration {} finished, {} particles remaining", iteration, particles.size());
            current_time += delta_t;
        }
        SPDLOG_INFO("Simulation completed: {} iterations, {} particles remaining", iteration, particles.size());
    }
};

}  // namespace mol_sim

#endif
