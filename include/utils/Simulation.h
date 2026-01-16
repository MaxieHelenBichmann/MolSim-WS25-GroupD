#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "exceptions/SimulationException.h"
#include "io/CheckpointWriter.h"
#include "io/OutputWriter.h"
#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
#include "particles/container/LinkedCellContainer.h"
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
     * @brief Force type, used in checkpointing.
     */
    Force force;
    /**
     * @brief Writer used for output.
     */
    const OutputWriter& writer;
    /**
     * @brief Writer used for output.
     */
    const CheckpointWriter& cp_writer;
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
    size_t frequency_output;

    /**
     * @brief Frequency of checkpoint writing.
     */
    size_t frequency_checkpoint;

    /**
     * @brief Base name for output files.
     */
    std::string base_name;

    /**
     * @brief Number of dimensions in the simulation.
     */
    size_t dimensions;
    /**
     * @brief Cutoff radius for particles in proximity.
     */
    double cutoff_radius;
    /**
     * @brief Total energy of the current system.
     */
    double total_energy;
    /**
     * @brief Target temperature of the system.
     */
    double target_temp;
    /**
     * @brief maximum allowed temperature change in one iteration.
     */
    double delta_temp;
    /**
     * @brief Frequency with which the thermostat is applied.
     */
    size_t thermostat_freq;
    /**
     * @brief Graviational constant to be used for simulating gravity in the simulation.
     * Will be used to apply a force m * g_grav (y-axis) to each particle.
     * Make sure it's 0 if you do NOT want to simulate gravitational pull along y-axis.
     * Make sure it's NEGATIVE if you want the particle to be pulled DOWN the y-axis.
     */
    double g_grav;
    /**
     * @brief Flag if thermostat is enabled for this simulation.
     */
    bool thermo;
    /**
     * @brief Temporary storage for new particles created during boundary application.
     */
    std::vector<Particle> new_particles;

   public:
    /**
     * @brief Construct a new Simulation object and prepare for run() call.
     *
     * @param particles Container of particles to be simulated.
     * @param force_source Force source to be used in the simulation.
     * @param settings Simulation parameters.
     */
    Simulation(containerType& particles, const ForceSource& force_source, SettingsParam& settings,
               const OutputWriter& writer, const CheckpointWriter& cp_writer)
        : domain(std::move(settings.domain)),
          particles(particles),
          force_source(force_source),
          force(settings.force),
          writer(writer),
          cp_writer(cp_writer),
          delta_t(settings.delta_t),
          start_time(settings.start_time),
          end_time(settings.end_time),
          frequency_output(settings.frequency_output),
          frequency_checkpoint(settings.frequency_checkpoint),
          base_name(settings.base_name),
          dimensions(settings.dimensions),
          cutoff_radius(settings.cutoff),
          target_temp(settings.target_temp),
          delta_temp(settings.delta_temp),
          thermostat_freq(settings.thermostat_freq),
          g_grav(settings.g_grav),
          thermo(settings.thermo) {
        #pragma omp parallel for
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = (*it);
            #pragma omp atomic
            total_energy += p.getM() * R3::scalarProduct(p.getV(), p.getV());
        }
        total_energy *= 0.5;
    }
/**
 * TODO:
 * maybe take into account the relation between creation and #iterations per thread
 * e.g. at 1400 particles and 50 bigger overhead for thread creation than thread
 * computing. In that case sweet spot at lower number of threads
 */

    double& getTotalEnergy() { return total_energy; }

    /**
     * @brief Removes all particles in the Halo from the container.
     */
    void removeParticles() {
        // Collect indices of particles to remove using halo iterator
        SPDLOG_DEBUG("Container has currently {} particles before erase", particles.size());
        std::vector<size_t> to_remove;
        #pragma omp parallel for 
        for (auto it = particles.haloBegin(); it != particles.haloEnd(); ++it) {
            size_t idx = &(*it) - &particles[0];
            to_remove.push_back(idx);
        }

        SPDLOG_DEBUG("Added {} (ghost) particles to remove", to_remove.size());

        // Sort in descending order to remove from end first (avoids index shifting issues that cause segfaults)
        std::sort(to_remove.begin(), to_remove.end(), std::greater<size_t>());  // NOLINT

        // Remove particles using the standard vector iterator version 
        #pragma omp parallel for
        for (size_t idx : to_remove) {
            particles.eraseParticle(particles.begin() + static_cast<std::ptrdiff_t>(idx));
        }
        SPDLOG_DEBUG("Container has currently {} particles after erase", particles.size());
    }

    /**
     * @brief Applies the necessary boundary conditions to the particles.
     */
    void applyBoundaries() {
//============================ISSUE (it = updatePartPos)======================================
        for (auto it = particles.begin(); it != particles.end();) {
            (*it).getOldF() = (*it).getF();
            (*it).getF() = Vector<double, 3>();
            // TODO: Optimization to only call this for relevant particles
            //std::vector<Particle> tmp = domain.applyBoundary(*it, force_source);
            //#pragma omp parallel for
            for (auto& p : domain.applyBoundary(*it, force_source)) {
                new_particles.push_back(p);
            }
            (*it).getMirrorLocations() = 0;
            // TODO: bit of an ugly workaround for now.
            R3 new_position = (*it).getX();
            (*it).getX() = (*it).getOldX();
            it = particles.updateParticlePosition(
                it, new_position);  // for now SimpleContainer + Periodic (and also Reflecting) needs this here
        }
//============================================================================================
    }

    /**
     * @brief Calculates the forces of every particle for the next time step.
     */
    void calculateF_LCC() {
        #pragma omp parallel for
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            Particle& p1 = *it;
            p1.getF()[1] += p1.getM() * g_grav;  // add gravitational pull along y-axis

            LinkedCellContainer::proximity_iterator<Particle, Cell> it_prox;
            LinkedCellContainer::proximity_iterator<Particle, Cell> it_prox_end;
            #pragma omp critical
            {
            it_prox = particles.proximityBegin(p1.getX(), it - particles.begin());
            it_prox_end = particles.proximityEnd(p1.getX());
            }

            //#pragma omp parallel for 
            //(disabled for now cause of parallelization overhead)
            //could be changed to #pragma omp parallel for if (cutoff > some_value) or smn like that
            for (; it_prox != it_prox_end; ++it_prox) {
                // Apply force directly (Newton's 3rd law: equal and opposite)
                #pragma omp critical
                {
                Particle& p2 = *it_prox;
                Vector<double, 3> force = force_source.applyForce(p1, p2);
                p1.getF() += force;
                p2.getF() -= force;
                }   
            }
        }
        // Calculate forces from mirrored/ghost particles
        #pragma omp parallel for
        for (const Particle& p1 : new_particles) {
            R3 lookup_pos = p1.getX();

            // For mirrored particles from periodic boundaries (type==1), their position is slightly
            // outside the domain. We need to clamp it to just inside the boundary region to find
            // the correct neighbors while preserving the actual position for force calculation.
            if (p1.getType() == 1) {  // Mirrored particle from periodic boundary
                R3 domain_size = domain.getDimension();
                constexpr double epsilon = 1e-6;  // Small offset to stay inside domain
                for (size_t dim = 0; dim < 3; ++dim) {
                    // Mirror slightly left of domain (x < 0) → clamp to just inside left boundary
                    if (lookup_pos[dim] < 0) {
                        lookup_pos[dim] = epsilon;
                    }
                    // Mirror slightly right of domain (x > domain) → clamp to just inside right boundary
                    else if (lookup_pos[dim] > domain_size[dim]) {
                        lookup_pos[dim] = domain_size[dim] - epsilon;
                    }
                }
            }

            auto it_prox = particles.proximityBegin(lookup_pos, particles.size());
            auto it_prox_end = particles.proximityEnd(lookup_pos);
            //#pragma omp parallel for 
            //(disabled for now cause of parallelization overhead)
            //could be changed to #pragma omp parallel for if (cutoff > some_value) or smn like that
            for (; it_prox != it_prox_end; ++it_prox) {
                Particle& p2 = *it_prox;
                Vector<double, 3> force = force_source.applyForce(p2, p1);
                p2.getF() = p2.getF() + force;
            }
        }
        new_particles.clear();
    }

    /**
     * @brief Calculates the forces of every particle for the next time step.
     */
    void calculateF() {
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            Particle& p1 = *it;
            p1.getF()[1] += p1.getM() * g_grav;  // add gravitational pull along y-axis

            auto it_prox = particles.proximityBegin(p1.getX(), it - particles.begin());
            auto it_prox_end = particles.proximityEnd(p1.getX());   

            //#pragma omp parallel for 
            //(disabled for now cause of parallelization overhead)
            //could be changed to #pragma omp parallel for if (cutoff > some_value) or smn like that
            for (; it_prox != it_prox_end; ++it_prox) {
                // Apply force directly (Newton's 3rd law: equal and opposite)
                Particle& p2 = *it_prox;
                Vector<double, 3> force = force_source.applyForce(p1, p2);
                p1.getF() += force;
                p2.getF() -= force;
            }   
        }
        // Calculate forces from mirrored/ghost particles
        for (const Particle& p1 : new_particles) {
            R3 lookup_pos = p1.getX();

            // For mirrored particles from periodic boundaries (type==1), their position is slightly
            // outside the domain. We need to clamp it to just inside the boundary region to find
            // the correct neighbors while preserving the actual position for force calculation.
            if (p1.getType() == 1) {  // Mirrored particle from periodic boundary
                R3 domain_size = domain.getDimension();
                constexpr double epsilon = 1e-6;  // Small offset to stay inside domain
                for (size_t dim = 0; dim < 3; ++dim) {
                    // Mirror slightly left of domain (x < 0) → clamp to just inside left boundary
                    if (lookup_pos[dim] < 0) {
                        lookup_pos[dim] = epsilon;
                    }
                    // Mirror slightly right of domain (x > domain) → clamp to just inside right boundary
                    else if (lookup_pos[dim] > domain_size[dim]) {
                        lookup_pos[dim] = domain_size[dim] - epsilon;
                    }
                }
            }

            auto it_prox = particles.proximityBegin(lookup_pos, particles.size());
            auto it_prox_end = particles.proximityEnd(lookup_pos);
            //#pragma omp parallel for 
            //(disabled for now cause of parallelization overhead)
            //could be changed to #pragma omp parallel for if (cutoff > some_value) or smn like that
            for (; it_prox != it_prox_end; ++it_prox) {
                Particle& p2 = *it_prox;
                Vector<double, 3> force = force_source.applyForce(p2, p1);
                p2.getF() = p2.getF() + force;
            }
        }
        new_particles.clear();
    }

    /**
     * @brief Calculates the positions of every particle for the next time step.
     */
    void calculateX() {
        #pragma omp parallel for
        for (auto& p : particles) {
            p.getOldX() = p.getX();
            p.getX() += (delta_t * p.getV()) + ((0.5 * delta_t * delta_t / p.getM()) * p.getF());
        }
    }

    /**
     * @brief Calculates the velocities of every particle for the next time step.
     */
    void calculateV(double scalar_factor) {
        double curr_energy = 0;
        #pragma omp parallel for
        for (auto& p : particles) {
            R3 new_v = scalar_factor * (p.getV() + ((0.5 * delta_t / p.getM()) * (p.getOldF() + p.getF())));
            p.getV() = new_v;
            curr_energy += p.getM() * R3::scalarProduct(new_v, new_v);
        }
        total_energy = 0.5 * curr_energy;
    }

    /**
     * @brief      Calculates the thermostat factor used to modulate velocity.
     *
     * @return     The thermostat factor.
     */
    double calculateThermostatFactor() {
        double curr_temp = (2.0 * total_energy) / (dimensions * particles.size());
        if (curr_temp == 0) {
            return 1;
        }
        double clamped_target = curr_temp + std::clamp((target_temp - curr_temp), -delta_temp, delta_temp);
        return sqrt(clamped_target / curr_temp);
    }

    /**
     * @brief Performs a full simulation run.
     * @throws SimulationException if an error occurs during output writing.
     */
    void run() {
        double current_time = start_time;
        [[maybe_unused]] int iteration = 0;

#ifdef ENABLE_CHECKPOINTING
        SettingsParam cp_settings;
        cp_settings.delta_t = delta_t;
        cp_settings.end_time = end_time;
        cp_settings.start_time = current_time;
        cp_settings.base_name = base_name;
        cp_settings.force = force;
        if constexpr (std::is_same_v<std::remove_cvref_t<containerType>, LinkedCellContainer>) {
            cp_settings.container_type = "LINKED";
        } else {
            cp_settings.container_type = "SIMPLE";
        }
        cp_settings.frequency_output = frequency_output;
        cp_settings.frequency_checkpoint = frequency_checkpoint;
        cp_settings.cutoff = cutoff_radius;
        cp_settings.target_temp = target_temp;
        cp_settings.delta_temp = delta_temp;
        cp_settings.thermostat_freq = thermostat_freq;
        cp_settings.dimensions = dimensions;
        cp_settings.g_grav = g_grav;
        cp_settings.thermo = thermo;
        auto cp_n = static_cast<size_t>(std::ceil((end_time - start_time) / delta_t));
#endif

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
            if constexpr (std::is_same_v<std::remove_cvref_t<containerType>, LinkedCellContainer>) {
                calculateF_LCC();
            } else {
                calculateF();
            }

            // 5. Calculate thermostat factor
            double thermo_factor = 1.0;
            // TODO: Iteration > 0 is a fix for the tests, discussion needed.
            if (thermo && iteration > 0 && iteration % thermostat_freq == 0) {
                thermo_factor = calculateThermostatFactor();
            }
            // 6. Calculate new velocities
            calculateV(thermo_factor);

            iteration++;
#ifdef ENABLE_IO
            if (iteration % frequency_output == 0) {
                try {
                    std::string out_name = base_name;
                    writer.plotParticles(particles, out_name, iteration);
                } catch (const std::runtime_error& e) {
                    SPDLOG_ERROR("Failed to plot particles at iteration {}: {}", iteration, e.what());
                    throw SimulationException("Error while plotting Particles: " + std::string(e.what()));
                }
            }
#endif
#ifdef ENABLE_CHECKPOINTING
            if (iteration % frequency_checkpoint == 0) {
                try {
                    cp_settings.start_time = current_time;
                    cp_writer.createCheckpoint(cp_settings, domain, particles, iteration, cp_n);
                } catch (const std::runtime_error& e) {
                    SPDLOG_ERROR("Failed to create a checkpoint at iteration {}: {}", iteration, e.what());
                    throw SimulationException("Error while creating checkpoint: " + std::string(e.what()));
                }
            }
#endif

            SPDLOG_DEBUG("Iteration {} finished, {} particles remaining", iteration, particles.size());
            current_time += delta_t;
        }
        SPDLOG_INFO("Simulation completed: {} iterations, {} particles remaining", iteration, particles.size());
#ifdef ENABLE_CHECKPOINTING
        try {
            cp_settings.start_time = current_time;
            cp_writer.createCheckpoint(cp_settings, domain, particles, iteration, cp_n);
        } catch (const std::runtime_error& e) {
            SPDLOG_ERROR("Failed to create a checkpoint at iteration {}: {}", iteration, e.what());
            throw SimulationException("Error while creating checkpoint: " + std::string(e.what()));
        }
#endif
    }
};

}  // namespace mol_sim

#endif
