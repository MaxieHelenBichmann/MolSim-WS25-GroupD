#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>

#include <memory>

#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/ParticleContainer.h"
#include "physics/ForceSource.h"
#include "physics/GravitationalForce.h"
#include "physics/LennardJonesForce.h"
#include "utils/Logging.h"
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
template <ParticleContainer containerType>
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
    std::unique_ptr<ForceSource> force_source;
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

    /**
     * @brief Calculates the forces of every particle for the next time step.
     * Calculates the forces of every particle. for the next time step. Using the specified force source and delta_t.
     */
    void calculateF() {
        std::vector<Vector<double, 3>> forces(particles.size(), Vector<double,3>());

        for (auto& p : particles) {
            p.getOldF() = p.getF();
            p.getF() = Vector<double, 3>();
        }

        for (size_t i = 0; i <= particles.size(); i++) {
            Particle p1 = particles[i];
            // compute row
            for (size_t j = i+1; j < particles.size(); j++) {
                Particle p2 = particles[j];
                Vector<double, 3> force = force_source->calculateForce(p1, p2);
                forces[j] = force;
                p1.getF() = p1.getF() + force;
            }
            // apply row in column 
            for (size_t j = i+1; j < particles.size(); j++) {
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
    void calculateFUnoptimized() {
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

   public:
    /**
     * @brief Construct a new Simulation object and prepare for run() call
     * This class implements a Builder Pattern, meaning all parameters need to be set before the run() call, which will
     * run the simulation.
     * @param particles Container of particles to be used in the simulation.
     * @param settings Simulation parameters. If relevant values are not set their default values in
     * include/utils/Default.h will be used instead.
     */
    Simulation(containerType& particles, SettingsParam& settings) : particles(particles) {
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        delta_t = settings.delta_t.value();
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        start_time = settings.start_time.value();
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        end_time = settings.end_time.value();
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        switch (settings.force_type.value()) {
            case GRAVITATIONAL:
                this->force_source = std::make_unique<GravitationalForce>();
                break;
            case LENNARDJONES:
                this->force_source = std::make_unique<LennardJonesForce>(
                    settings.epsilon.value(), settings.sigma.value());  // NOLINT(bugprone-unchecked-optional-access)
            default:
                break;
        }
    }

    /**
     * HOT FIX FOR Simulation_test.cpp. NEED TO CHANGE THAT FILE CAUSE WE CANT KEEP THIS CONSTRUCTOR (its ugly)
     */
    Simulation(containerType& particles, Force forceType, double delta_t, double start_time, double end_time)
        : particles(particles), delta_t(delta_t), start_time(start_time), end_time(end_time) {
        switch (forceType) {
            case GRAVITATIONAL:
                this->force_source = std::make_unique<GravitationalForce>();
                break;
            default:
                break;
        }
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
            if (iteration % 10 == 0) {
                try {
#ifdef ENABLE_VTK_OUTPUT
                    std::string out_name("MD_vtk");
                    VTKWriter writer;
#else
                    std::string out_name("MD_xyz");
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
