#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>

#include <cstddef>
#include <limits>

#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
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
 */
template <ParticleContainer containerType, ForceSource forceType>
class Simulation {
   private:
    /**
     * @brief Container of Particles to simulate.
     * Container of Particles to simulate. Container type is templated to work with our container concept.
     */
    containerType& particles;
    /**
     * @brief Pointer to our force source.
     * Pointer to our force source, for easy switching, force Source determined by forceType in constructor.
     */
    forceType& force_source;
    /**
     * @brief Time step of simulation.
     * Default value is 0.014.
     */
    double delta_t;
    /**
     * @brief Start time of simulation.
     * Default value is 0.
     */
    double start_time;
    /**
     * @brief End time of simulation.
     * Default value is 1000.
     */
    double end_time;

    size_t frequency;
    std::string base_name;

    /**
     * @brief Cutoff radius for particles in proximity.
     * Default value is infinity.
     */
    double cutoff_radius = std::numeric_limits<double>::infinity();

   public:
    /**
     * @brief Calculates the forces of every particle for the next time step.
     * Calculates the forces of every particle. for the next time step. Using the specified force source and delta_t.
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
     * Calculates the forces of every particle for the next time step, specified by delta_t.
     */
    void calculateX() {
        for (auto& p : particles) {
            p.getX() = p.getX() + (delta_t * p.getV()) + ((0.5 * delta_t * delta_t / p.getM()) * p.getF());
        }
    }

    /**
     * @brief Calculates the velocities of every particle for the next time step.
     * Calculates the forces of every particle for the next time step, specified by delta_t.
     */
    void calculateV() {
        for (auto& p : particles) {
            p.getV() = p.getV() + ((0.5 * delta_t / p.getM()) * (p.getOldF() + p.getF()));
        }
    }

    /**
     * @brief Construct a new Simulation object and prepare for run() call
     * This class implements a Builder Pattern, meaning all parameters need to be set before the run() call, which will
     * run the simulation.
     * @param particles Container of particles to be used in the simulation.
     * @param force_source Force source to be used in the simulation.
     * @param settings Simulation parameters. If relevant values are not set their default values in
     * include/utils/Default.h will be used instead.
     */
    Simulation(containerType& particles, forceType& force_source, SettingsParam& settings)
        : particles(particles), force_source(force_source) {
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        delta_t = settings.delta_t.value();
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        start_time = settings.start_time.value();
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        end_time = settings.end_time.value();
        frequency = settings.frequency.value();
        base_name = settings.base_name.value();
    }

    /**
     * @brief Performs a full simulation run.
     * Performs a full simulation run, using the specified delta_t and end_time.
     */
    void run() {
        double current_time = start_time;
        [[maybe_unused]] int iteration = 0;

        // for this loop, we assume: current x, current f and current v are known
        while (current_time < end_time) {
            // calculate new x
            calculateX();
            // calculate new f
            calculateF();
            // calculate new v
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
