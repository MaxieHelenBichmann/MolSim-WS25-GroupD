#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>
#include <memory>

#include "exceptions/SimulationException.h"
#include "io/CheckpointWriter.h"
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
    bool thermo;
    /**
     * @brief Vector of particle indicies that should be removed the next time we call 'removeParticles'
     */
    std::vector<size_t> to_remove;

    std::unique_ptr<LinkedCellContainer> checker;

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
        for (const Particle& p : particles) {
            total_energy += p.getM() * R3::scalarProduct(p.getV(), p.getV());
        }
        total_energy *= 0.5;
        checker = std::make_unique<LinkedCellContainer>(domain.getDimension(), cutoff_radius);
    }

    double& getTotalEnergy() { return total_energy; }

    /**
     * @brief Calculates the positions of every particle for the next time step.
     */
    void calculateX() {
        for (auto& p : particles) {
            if (p.getType() == 1) { // mirror particles will be handled in applyBoundaries
                p.getOldX() = p.getX();
                Particle parent = particles[p.getParent()];
                R3 delta_x = parent.getX() - parent.getOldX();
                p.getX() = p.getX() + delta_x;
            } else {
                p.getOldX() = p.getX();
                p.getX() = p.getX() + (delta_t * p.getV()) + ((0.5 * delta_t * delta_t / p.getM()) * p.getF());
            }
        }
    }

    bool isInHalo(Particle& p) {
        //TODO: Optimization possible. If p is exactly on boundary this function
        //will return false, even though it should return true. Since this function (so far)
        //is only used for updating mirror particles, this will lead to mirror particles that are exactly
        //on a boundary being unnecessarily deleted and recreated.
        return checker->fitsContainer(p.getX()) && !checker->fitsDomain(p.getX());
    }

    void updateMirrorParticles(Particle& p, const R3& delta_x) { 
        auto begin = p.getMirrorParticles().begin();
        auto end = p.getMirrorParticles().end();
        for (auto mirror_idx = begin; mirror_idx != end;) {
            Particle mirrorParticle = particles[*mirror_idx];
            mirrorParticle.getX() = mirrorParticle.getX() + delta_x;
            if (!isInHalo(mirrorParticle)) {
                to_remove.push_back(*mirror_idx);
                mirror_idx = p.getMirrorParticles().erase(mirror_idx);
            } else {
                ++mirror_idx;
            }
        }
    }

    /**
     * @brief Applies the necessary boundary conditions to the particles.
     */
    void applyBoundaries() {
        for (auto it = particles.begin(); it != particles.end();) {
            if ((*it).getType() == 1) { // don't do any of this for mirror particles
                break; //fine assuming after we read first mirror particle there are only other mirror particles that follow
            }
            (*it).getOldF() = (*it).getF();
            (*it).getF() = Vector<double, 3>();
            // TODO: Optimization to only call this for relevant particles
            std::vector<Particle> new_particles;
            R3 positionBeforeBoundary = (*it).getX();
            for (auto& p : domain.applyBoundary(*it, force_source)) {
                new_particles.push_back(p);
                if (p.getType() == 1) {
                    p.getParent() = (&(*it) - &particles[0]);
                }
            }
            if (!checker->fitsContainer((*it).getX())) { // if new position OOB remove this particle and its mirrors
                to_remove.push_back(&(*it) - &particles[0]); 
                auto begin = (*it).getMirrorParticles().begin();
                auto end = (*it).getMirrorParticles().end();
                for (auto mirror_idx = begin; mirror_idx != end; ++mirror_idx) {
                    to_remove.push_back(*mirror_idx);
                }
                continue;
            } else {
                //1) update positions of already existing mirror particles
                updateMirrorParticles((*it), (*it).getX() - positionBeforeBoundary);
                //2) add the newly generated mirror particles (they're already in the correct position)
                for (const auto& p : new_particles) {
                    particles.addParticle(p);
                    particles[p.getParent()].getMirrorParticles().push_back(particles.size()-1);
                }
            }
            // TODO: bit of an ugly workaround for now.
            R3 new_position = (*it).getX();
            (*it).getX() = (*it).getOldX();
            it = particles.updateParticlePosition(
                it, new_position);  // for now SimpleContainer + Periodic (and also Reflecting) needs this here
        }
    }

    /**
     * @brief Removes all particles in the Halo from the container.
     */
    void removeParticles() {
        // Collect indices of particles to remove using halo iterator
        SPDLOG_DEBUG("Container has currently {} particles before erase", particles.size());
        for (auto it = particles.haloBegin(); it != particles.haloEnd(); ++it) {
            //0) don't remove mirror particles if parent isn't also going to be removed here.
            //lone mirror particles that should be removed because they are no longer in the halo, 
            //are added to to_remove in calculateX
            if ((*it).getType() == 1) {
                continue;
            }
            //1) remove normal particle
            size_t idx = &(*it) - &particles[0];
            to_remove.push_back(idx);
            //2) remove its mirror particles
            auto begin = (*it).getMirrorParticles().begin();
            auto end = (*it).getMirrorParticles().end();
            for (auto mirror_idx = begin; mirror_idx != end; ++mirror_idx) {
                to_remove.push_back(*mirror_idx);
            }
        }

        SPDLOG_DEBUG("Added {} (ghost) particles to remove", to_remove.size());

        // Sort in ascending order so shifting mirror particle indicies works as intended
        std::sort(to_remove.begin(), to_remove.end(), std::less<size_t>());  // NOLINT
        
        // Shift mirror particle indicies
        // The following code assumes:
        // 1) that to_remove contains the indicies of *all* particles we'd like to remove
        //    including mirror particles, sorted in ascending order.
        // 2) that the normal particles (of type 0) are all at the start of the particles collection, 
        //    and the mirror particles (particles of type 1) are then all at the end of the particles collection.
        //    Meaning once we read the first type 0 particle, no particles of another type should follow until we read a particle
        //    of type 1 for the first time, in which case, no particles of a type other than 1 should follow.
        size_t shift = 0;
        size_t to_remove_idx = 0;
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            size_t idx = &(*it) - &particles[0];
            if (to_remove_idx < to_remove.size() && idx == to_remove[to_remove_idx]) {
                shift++;
                to_remove_idx++;
                continue;
            }
            if (shift > 0) {
                if ((*it).getType() == 0) { //tell mirror particles how parent shifted
                    auto begin = (*it).getMirrorParticles().begin();
                    auto end = (*it).getMirrorParticles().end();
                    for (auto mirror_it = begin; mirror_it != end; ++mirror_it) {
                        particles[*mirror_it].getParent() -= shift;
                    }
                } else if ((*it).getType() == 1) { //tell parent how mirror particles shifted
                    particles[(*it).getParent()].getMirrorParticles()[(*it).getMirrorIdx()] -= shift;
                }
            }
        }

        // Sort in descending order to remove from end first (avoids index shifting issues)
        std::sort(to_remove.begin(), to_remove.end(), std::greater<size_t>());  // NOLINT

        // Remove particles using the standard vector iterator version
        for (size_t idx : to_remove) {
            particles.eraseParticle(particles.begin() + static_cast<std::ptrdiff_t>(idx));
        }

        to_remove.clear();
        SPDLOG_DEBUG("Container has currently {} particles after erase", particles.size());
    }

    /**
     * @brief Calculates the forces of every particle for the next time step.
     */
    void calculateF() {
        size_t idx = 0;
        for (auto it = particles.begin(); it != particles.end(); ++it, idx++) {
            Particle& p1 = *it;
            if (p1.getType() == 1) {
                break; //fine assuming after we read first mirror particle there are only other mirror particles that follow
            }
            p1.getF()[1] += p1.getM() * g_grav;  // add gravitational pull along y-axis
            auto it_prox = particles.proximityBegin(p1.getX(), idx);
            auto it_prox_end = particles.proximityEnd(p1.getX());
            for (; it_prox != it_prox_end; ++it_prox) {
                Particle& p2 = *it_prox;
                Vector<double, 3> force = force_source.applyForce(p1, p2);
                // Apply force directly (Newton's 3rd law: equal and opposite)
                p1.getF() = p1.getF() + force;
                if (p2.getType() == 1) {
                    continue; //don't touch mirror particles
                }
                p2.getF() = p2.getF() - force;
            }
        }
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
                    cp_writer.createCheckpoint(domain, particles, iteration, force, delta_t, current_time, end_time,
                                               frequency_output, frequency_checkpoint, base_name, cutoff_radius,
                                               target_temp, delta_temp, thermostat_freq,
                                               static_cast<size_t>(std::ceil((end_time - start_time) / delta_t)));
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
    }
};

}  // namespace mol_sim

#endif
