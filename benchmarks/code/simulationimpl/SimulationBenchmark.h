#ifndef SIMULATION_BENCHMARK_H
#define SIMULATION_BENCHMARK_H

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include "io/OutputWriter.h"
#include "particles/ParticleContainer.h"
#include "particles/container/domain/Domain.h"
#include "physics/ForceSource.h"
#include "utils/Settings.h"

namespace mol_sim {

/**
 * @brief Class for running Simulation for the LinkedCellContainerDirect.
 * We need a specialized version to be able to fairly benchmark the different container implementations.
 * Implements all methods of the Simulation (so documentation is analogous), but only used for benchmarking.
 * It updates the particle positions in a different way and does not use the N3L optimization for force calculations.
 *
 * Only used for benchmarking SimpleContainerControl, LinkedCellContainerDirect and LinkedCellContainerExplicit.
 */
template <ParticleContainer containerType>
class SimulationBenchmark {
   private:
    Domain domain;
    containerType& particles;
    const ForceSource& force_source;
    const OutputWriter& writer;
    double delta_t;
    double start_time;
    double end_time;
    size_t frequency;
    std::string base_name;
    double cutoff_radius;

   public:
    SimulationBenchmark(containerType& particles, const ForceSource& force_source, SettingsParam& settings,
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
          cutoff_radius(settings.cutoff) {}

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

    void applyBoundaries() {
        for (auto it = particles.begin(); it != particles.end();) {
            (*it).getOldF() = (*it).getF();
            (*it).getF() = Vector<double, 3>();
            domain.applyBoundary(*it, force_source);
        }
    }

    void calculateF() {
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            for (auto it_prox = particles.proximityBegin((*it).getX(), 0);
                 it_prox != particles.proximityEnd((*it).getX()); ++it_prox) {
                (*it).getF() = (*it).getF() + force_source.applyForce(*it, *it_prox);
            }
        }
    }

    void calculateX() {
        for (auto it = particles.begin(); it != particles.end();) {
            R3 new_x =
                (*it).getX() + (delta_t * (*it).getV()) + ((0.5 * delta_t * delta_t / (*it).getM()) * (*it).getF());
            it = particles.updateParticlePosition(it, new_x);
        }
    }

    void calculateV() {
        for (auto& p : particles) {
            p.getV() = p.getV() + ((0.5 * delta_t / p.getM()) * (p.getOldF() + p.getF()));
        }
    }

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

            // 6. Calculate new velocities
            calculateV();

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
