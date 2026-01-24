#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>

#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "io/CheckpointWriter.h"
#include "io/OutputWriter.h"
#include "io/StatsWriter.h"
#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
#include "particles/container/domain/Domain.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "physics/singleforces/HarmonicForce.h"
#include "physics/singleforces/SingleForceSource.h"
#include "physics/targettedforces/TargetForceSource.h"
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
     * @brief Force sources for calculating particle interactions.
     */
    const std::vector<std::unique_ptr<PairwiseForceSource>>& pairwise_force_sources;

    const std::vector<std::unique_ptr<SingleForceSource>>& single_force_sources;

    std::vector<PairwiseForce> pairwise_forces;
    std::vector<SingleForce> single_forces;
    /**
     * @brief Writer used for output.
     */
    const OutputWriter& writer;
    /**
     * @brief Writer used for output.
     */
    const CheckpointWriter& cp_writer;
    /**
     * @brief Writer used for statistics.
     */
    const StatsWriter& stats_writer;

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
     * @brief Frequency of statistics writing.
     */
    size_t frequency_stats_diff;
    /**
     * @brief Frequency of statistics writing.
     */
    size_t frequency_stats_rdf;

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
     * @brief Flag if thermostat is enabled for this simulation.
     */
    bool thermo;
    /**
     * @brief Temporary storage for new particles created during boundary application.
     */
    std::vector<Particle> new_particles;

    bool target_force_enabled = false;

    TargetForceSource target_force;

    /**
     * @brief Gravitational acceleration vector for single GRAV force.
     */
    R3 g_grav_vec;

    /**
     * @brief Spring constant for HARMONIC force.
     */
    double k;

    /**
     * @brief Equilibrium distance for HARMONIC force.
     */
    double r_0;

    /**
     * @brief Parallelization strategy for force calculation.
     */
    ParallelizationStrategy strategy;

    /**
     * @brief Function pointer to the selected force calculation method.
     * Set once during construction to avoid runtime checks every iteration.
     */
    std::function<void(const size_t)> calculate_forces;

   public:
    /**
     * @brief Construct a new Simulation object and prepare for run() call.
     *
     * @param particles Container of particles to be simulated.
     * @param force_source Force source to be used in the simulation.
     * @param settings Simulation parameters.
     */
    Simulation(containerType& particles,
               const std::vector<std::unique_ptr<PairwiseForceSource>>& pairwise_force_sources,
               const std::vector<std::unique_ptr<SingleForceSource>>& single_force_sources, SettingsParam& settings,
               const OutputWriter& writer, const CheckpointWriter& cp_writer, const StatsWriter& stats_writer)
        : domain(std::move(settings.domain)),
          particles(particles),
          pairwise_force_sources(pairwise_force_sources),
          single_force_sources(single_force_sources),
          pairwise_forces(settings.pairwise_forces),
          single_forces(settings.single_forces),
          writer(writer),
          cp_writer(cp_writer),
          stats_writer(stats_writer),
          delta_t(settings.delta_t),
          start_time(settings.start_time),
          end_time(settings.end_time),
          frequency_output(settings.frequency_output),
          frequency_checkpoint(settings.frequency_checkpoint),
          frequency_stats_diff(settings.stats_freq_diffusion),
          frequency_stats_rdf(settings.stats_freq_rdf),
          base_name(settings.base_name),
          dimensions(settings.dimensions),
          cutoff_radius(settings.cutoff),
          target_temp(settings.target_temp),
          delta_temp(settings.delta_temp),
          thermostat_freq(settings.thermostat_freq),
          thermo(settings.thermo),
          target_force_enabled(settings.target_force_enabled),
          target_force(settings.target_force_direction, settings.target_force_magnitude,
                       settings.target_force_max_iterations),
          g_grav_vec(settings.g_grav_vec),
          k(settings.k),
          r_0(settings.r_0),
          strategy(settings.strategy) {
#ifdef _OPENMP
        // Set function pointer based on strategy to avoid runtime checks every iteration
        if (strategy == ParallelizationStrategy::COLORING) {
            calculate_forces = [this](const size_t iteration) { this->calculateFColored(iteration); };
        } else {
            calculate_forces = [this](const size_t iteration) { this->calculateF(iteration); };
        }
#else
        calculate_forces = [this](const size_t iteration) { this->calculateF(iteration); };
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : total_energy) default(none) shared(particles)
#endif
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = (*it);
            // NOLINTNEXTLINE
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
        for (auto it = particles.haloBegin(); it != particles.haloEnd(); ++it) {
            size_t idx = &(*it) - &particles[0];
            to_remove.push_back(idx);
        }

        SPDLOG_DEBUG("Added {} (ghost) particles to remove", to_remove.size());

        // Sort in descending order to remove from end first (avoids index shifting issues that cause segfaults)
        std::sort(to_remove.begin(), to_remove.end(), std::greater<size_t>());  // NOLINT

        for (size_t idx : to_remove) {
            particles.eraseParticle(particles.begin() + static_cast<std::ptrdiff_t>(idx));
        }
        SPDLOG_DEBUG("Container has currently {} particles after erase", particles.size());
    }

    /**
     * @brief Applies the necessary boundary conditions to the particles.
     */
    void applyBoundaries() {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 100) default(none) shared(particles, domain)
#endif
        for (auto& p : particles) {
            p.getOldF() = p.getF();
            p.getF() = Vector<double, 3>();
            // TODO: Optimization to only call this for relevant particles
            domain.applyBoundary(p);
            p.getMirrorLocations() = 0;
        }

        for (auto it = particles.begin(); it != particles.end();) {
            R3 new_position = (*it).getX();
            (*it).getX() = (*it).getOldX();
            it = particles.updateParticlePosition(it, new_position);
            // for now SimpleContainer + Periodic (and also Reflecting) needs this here
        }
    }
    /**
     * @brief Calculates the forces of every particle for the next time step.
     */
    void calculateF(const size_t iteration) {
        particles.prepareForParallelIteration();

        const size_t num_particles = particles.size();
        const R3& domain_size = domain.getDimension();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 100) default(none)                                                          \
    shared(particles, single_force_sources, target_force, pairwise_force_sources, num_particles, target_force_enabled, \
               iteration, domain_size)
#endif

        for (size_t i = 0; i < num_particles; ++i) {
            Particle& p1 = particles[i];
            Vector<double, 3> f1_accumulated{};

            // 1. Single Forces
            for (const auto& force_source : single_force_sources) {
                f1_accumulated += force_source->applyForce(p1);
            }

            // 2. Target Forces
            if (target_force_enabled && iteration < target_force.getMaxIterations()) {
                f1_accumulated += target_force.applyForce(p1);
            }

            // 3. Normal Pairwise Forces (Push to p2, Add to p1)
            const R3& p1_pos = p1.getX();
            auto it_prox = particles.proximityBegin(p1_pos, i);
            auto it_prox_end = particles.proximityEnd(p1_pos);

            for (; it_prox != it_prox_end; ++it_prox) {
                Particle& p2 = *it_prox;
                for (const auto& force_source : pairwise_force_sources) {
                    const Vector<double, 3> force = force_source->applyForce(p1, p2);
                    f1_accumulated += force;
                    p2.getF().atomicSubtract(force);
                }
            }

            // 4. Mirror Pairwise Forces (Gather to p1)
            for (const R3& mirr_pos : p1.getMirrorPositions()) {
                // Create a local ghost particle to avoid modifying p1's position in shared memory
                Particle p1_ghost = p1;
                p1_ghost.getX() = mirr_pos;
                R3 lookup_pos = mirr_pos;

                constexpr double epsilon = 1e-6;
                for (size_t dim = 0; dim < 3; ++dim) {
                    if (lookup_pos[dim] < 0) {
                        lookup_pos[dim] = epsilon;
                    } else if (lookup_pos[dim] > domain_size[dim]) {
                        lookup_pos[dim] = domain_size[dim] - epsilon;
                    }
                }
                auto it_prox_mirror = particles.proximityBegin(lookup_pos, particles.size());
                auto it_prox_mirror_end = particles.proximityEnd(lookup_pos);

                for (; it_prox_mirror != it_prox_mirror_end; ++it_prox_mirror) {
                    Particle& p2 = *it_prox_mirror;
                    for (const auto& force_source : pairwise_force_sources) {
                        const Vector<double, 3> force = force_source->applyForce(p2, p1_ghost);
                        // Newton's 3rd Law: Force on p1 is -Force on p2.
                        // We gather the force on p1 instead of scattering to p2 to avoid atomic ops on p2.
                        f1_accumulated -= force;
                    }
                }
            }

            // Final Apply to p1
            p1.getF().atomicAdd(f1_accumulated);
        }
    }
#ifdef _OPENMP
    /**
     * @brief Calculates forces using cell coloring to avoid atomic operations.
     *
     * Uses 8-color 3D checkerboard pattern where cells of same color never interact.
     * This eliminates all atomic operations, enabling perfect parallel scaling.
     * Only works with LinkedCellContainer.
     *
     * @param iteration Current simulation iteration number
     */
    void calculateFColored(const size_t iteration) {
        if constexpr (!std::is_same_v<containerType, LinkedCellContainer>) {
            SPDLOG_WARN("calculateFColored only works with LinkedCellContainer, falling back to calculateF");
            calculateF(iteration);
        } else {
            particles.prepareForParallelIteration();
            const R3& domain_size = domain.getDimension();

            // Get cell grid dimensions from LinkedCellContainer
            auto& lcc = static_cast<LinkedCellContainer&>(particles);
            const auto& num_cells = lcc.getNumCells();

#pragma omp parallel default(none)                                                                            \
    shared(particles, single_force_sources, target_force, pairwise_force_sources, dimensions, lcc, num_cells, \
               target_force_enabled, iteration, domain_size)
            {
                // 8-color 3D checkerboard: process each color sequentially, parallelize within color
                size_t num_colors = dimensions == 3 ? 8 : 4;
                for (size_t color = 0; color < num_colors; ++color) {
                    const size_t color_i = color & 1;         // bit 0
                    const size_t color_j = (color >> 1) & 1;  // bit 1
                    const size_t color_k = (color >> 2) & 1;  // bit 2

#pragma omp for schedule(dynamic) collapse(3)
                    for (size_t ci = color_i; ci < num_cells[0]; ci += 2) {
                        for (size_t cj = color_j; cj < num_cells[1]; cj += 2) {
                            for (size_t ck = color_k; ck < num_cells[2]; ck += 2) {
                                const size_t cell_idx = lcc.cellIndex(ci, cj, ck);
                                const Cell& cell = lcc.getCell(cell_idx);

                                // Process all particles in this cell
                                auto it = cell.stableIteratorBegin();
                                auto it_end = cell.stableIteratorEnd();
                                for (; it != it_end; ++it) {
                                    const size_t p1_idx = *it;
                                    Particle& p1 = particles[p1_idx];

                                    // 1. Single Forces
                                    for (const auto& force_source : single_force_sources) {
                                        p1.getF() += force_source->applyForce(p1);
                                    }

                                    // 2. Target Forces
                                    if (target_force_enabled && iteration < target_force.getMaxIterations()) {
                                        p1.getF() += target_force.applyForce(p1);
                                    }

                                    // 3. Normal Pairwise Forces
                                    const R3& p1_pos = p1.getX();
                                    auto it_prox = particles.proximityBegin(p1_pos, p1_idx);
                                    auto it_prox_end = particles.proximityEnd(p1_pos);

                                    for (; it_prox != it_prox_end; ++it_prox) {
                                        Particle& p2 = *it_prox;
                                        for (const auto& force_source : pairwise_force_sources) {
                                            const Vector<double, 3> force = force_source->applyForce(p1, p2);
                                            p1.getF() += force;
                                            p2.getF().atomicSubtract(force);
                                        }
                                    }

                                    // 4. Mirror Pairwise Forces
                                    for (const R3& mirr_pos : p1.getMirrorPositions()) {
                                        Particle p1_ghost = p1;
                                        p1_ghost.getX() = mirr_pos;
                                        R3 lookup_pos = mirr_pos;

                                        constexpr double epsilon = 1e-6;
                                        for (size_t dim = 0; dim < 3; ++dim) {
                                            if (lookup_pos[dim] < 0) {
                                                lookup_pos[dim] = epsilon;
                                            } else if (lookup_pos[dim] > domain_size[dim]) {
                                                lookup_pos[dim] = domain_size[dim] - epsilon;
                                            }
                                        }
                                        auto it_prox_mirror = particles.proximityBegin(lookup_pos, particles.size());
                                        auto it_prox_mirror_end = particles.proximityEnd(lookup_pos);

                                        for (; it_prox_mirror != it_prox_mirror_end; ++it_prox_mirror) {
                                            Particle& p2 = *it_prox_mirror;
                                            for (const auto& force_source : pairwise_force_sources) {
                                                const Vector<double, 3> force = force_source->applyForce(p1_ghost, p2);
                                                p1.getF() += force;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif
    /**
     * @brief Calculates the positions of every particle for the next time step.
     */
    void calculateX() {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(particles, delta_t)
#endif
        for (auto& p : particles) {
            p.getMirrorPositions().clear();
            p.getOldX() = p.getX();
            p.getX() += (delta_t * p.getV()) + ((0.5 * delta_t * delta_t / p.getM()) * p.getF());
        }
    }

    /**
     * @brief Calculates the velocities of every particle for the next time step.
     */
    void calculateV(const double scalar_factor) {
        double curr_energy = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : curr_energy) default(none) \
    shared(particles, scalar_factor, delta_t)
#endif
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
        for (const auto& source : single_force_sources) {
            if (source->getType() == SingleForce::HARMONIC) {
                static_cast<HarmonicForce*>(source.get())->setContainer(particles);
            }
        }

#ifdef ENABLE_CHECKPOINTING
        SettingsParam cp_settings;
        cp_settings.delta_t = delta_t;
        cp_settings.end_time = end_time;
        cp_settings.start_time = current_time;
        cp_settings.base_name = base_name;
        cp_settings.pairwise_forces = pairwise_forces;
        cp_settings.single_forces = single_forces;
        cp_settings.g_grav_vec = g_grav_vec;
        cp_settings.k = k;
        cp_settings.r_0 = r_0;
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
        cp_settings.thermo = thermo;
        cp_settings.strategy = strategy;
        cp_settings.target_force_enabled = target_force_enabled;
        if (target_force_enabled) {
            cp_settings.target_force_direction = target_force.getDirection();
            cp_settings.target_force_magnitude = target_force.getMagnitude();
            cp_settings.target_force_max_iterations = target_force.getMaxIterations();
        }
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
            calculate_forces(iteration);

            // 5. Calculate thermostat factor
            double thermo_factor = 1.0;
            if (thermo && iteration % thermostat_freq == 0) {
                thermo_factor = calculateThermostatFactor();
            }
            // 6. Calculate new velocities
            calculateV(thermo_factor);

            iteration++;
#ifdef ENABLE_IO
            // 7. Write output (optional)
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
            // 8. Write checkpoint (optional)
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
#ifdef ENABLE_STATS
            // 9. Write diffusion statistics (optional)
            if (iteration % frequency_stats_diff == 0) {
                try {
                    stats_writer.plotDiffusion(particles, iteration);
                } catch (const std::runtime_error& e) {
                    SPDLOG_ERROR("Failed to plot diffusion at iteration {}: {}", iteration, e.what());
                }
            }
            // 10. Write RDF statistics (optional)
            if (iteration % frequency_stats_rdf == 0) {
                try {
                    stats_writer.plotRDF(particles, iteration);
                } catch (const std::runtime_error& e) {
                    SPDLOG_ERROR("Failed to plot RDF at iteration {}: {}", iteration, e.what());
                }
            }
#endif

            SPDLOG_DEBUG("Iteration {} finished, {} particles remaining", iteration, particles.size());
            current_time += delta_t;
        }
        SPDLOG_INFO("Simulation completed: {} iterations, {} particles remaining", iteration, particles.size());
#ifdef ENABLE_CHECKPOINTING
        // Final checkpoint at end of simulation (optional)
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
