#ifndef SIMULATION_H
#define SIMULATION_H

#include <spdlog/spdlog.h>

#include <cstddef>
#include <limits>
#include <memory>
#include <typeinfo>

#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
#include "particles/boundaries/BoundaryCondition.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Reflecting.h"
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
 */
template <ParticleContainer containerType, ForceSource forceType>
class Simulation {
   private:
    /**
     * @brief Simulation domain.
     * 
     */
    Domain<containerType>& domain;
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
     * @deprecated This constructor only exists as a hot fix for the benchmarks 
     * after my (Georg) refactorings to the codebase.
     */
    Simulation(containerType& particles, forceType& force_source, SettingsParam& settings)
        : particles(particles), force_source(force_source), 
          domain(std::get<Domain<containerType>>(settings.domain)) {
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
     * @brief Construct a new Simulation object and prepare for run() call
     * This class implements a Builder Pattern, meaning all parameters need to be set before the run() call, which will
     * run the simulation.
     * @param settings Simulation parameters. If relevant values are not set their default values in
     * include/utils/Settings.h will be used instead.
     * @param force_source Force source to be used in the simulation.
     */
    Simulation(SettingsParam& settings, forceType& force_source)
        : force_source(force_source),
          domain(std::get<Domain<containerType>>(settings.domain)) {
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        delta_t = settings.delta_t.value();
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        start_time = settings.start_time.value();
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        end_time = settings.end_time.value();
        frequency = settings.frequency.value();
        base_name = settings.base_name.value();
        particles = domain.getParticles(); 
    }

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
     * Calculates the forces of every particle for the next time step, specified by delta_t, for the provided
     * container.
     */
    void calculateX() {
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            const auto new_position =
                (*it).getX() + (delta_t * (*it).getV()) + ((0.5 * delta_t * delta_t / (*it).getM()) * (*it).getF());
            particles.updateParticlePosition(it, new_position);
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

    void applyBoundary() {
        domain.getBoundary(BoundaryLocation::LEFT)->applyBoundary();  
        domain.getBoundary(BoundaryLocation::RIGHT)->applyBoundary();  
        domain.getBoundary(BoundaryLocation::UPPER)->applyBoundary();  
        domain.getBoundary(BoundaryLocation::LOWER)->applyBoundary();  
        domain.getBoundary(BoundaryLocation::FRONT)->applyBoundary();  
        domain.getBoundary(BoundaryLocation::BACK)->applyBoundary();  
    }

    /**
     * TODO: Implement analogy to halo cell and boundary cell iterator in SimpleContainer.
     * (maybe even add those to the concept) then uncomment the code inside this function.
     * 
     * @brief Clears the halo cells, i.e. removes all particles that are beyond the specified simulation domain.
     */
    void cleanBoundary() {
/*         for (auto& it = particles.haloBegin(); it != particles.haloEnd(); it++) {
            particles.eraseParticle(it);
        } */
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
            cleanBoundary();
            applyBoundary();
            calculateF();
            cleanBoundary();
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
