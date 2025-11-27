#include <spdlog/spdlog.h>

#include "io/CLIParse.h"
#include "io/fileReader/YAMLReaderException.h"
#include "particles/container/SimpleContainer.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/ContainerRef.h"
#include "physics/GravitationalForce.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

using namespace mol_sim;

void startSimulation(ContainerRef particles, SettingsParam& settings) {
    CONTAINER_REF particles_var = particles.getInstance();  
    if (std::holds_alternative<SimpleContainer*>(particles_var)) {
    switch (settings.force.value()) {
        case GRAVITATIONAL: {
            GravitationalForce grav_force;
            SimpleContainer deref_particles = *(std::get<SimpleContainer*>(particles_var)); 
            Simulation<SimpleContainer, GravitationalForce> simulation(deref_particles, settings, grav_force);
            simulation.run();
            return;
        }
        case LENNARDJONES: {
            LennardJonesForce lj_force;
            SimpleContainer deref_particles = *(std::get<SimpleContainer*>(particles_var)); 
            Simulation<SimpleContainer, LennardJonesForce> simulation(deref_particles, settings, lj_force);
            simulation.run();
            return;
        }
        default:
            SPDLOG_ERROR("Unrecognized force source!");
            return;
    }
    } else if (std::holds_alternative<LinkedCellContainer*>(particles_var)) {
    switch (settings.force.value()) {
        case GRAVITATIONAL: {
            GravitationalForce grav_force;
            LinkedCellContainer deref_particles = *(std::get<LinkedCellContainer*>(particles_var)); 
            Simulation<LinkedCellContainer, GravitationalForce> simulation(deref_particles, settings, grav_force);
            simulation.run();
            return;
        }
        case LENNARDJONES: {
            LennardJonesForce lj_force;
            LinkedCellContainer deref_particles = *(std::get<LinkedCellContainer*>(particles_var)); 
            Simulation<LinkedCellContainer, LennardJonesForce> simulation(deref_particles, settings, lj_force);
            simulation.run();
            return;
        }
        default:
            SPDLOG_ERROR("Unrecognized force source!");
            return;
        }
    } 
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char* argsv[]) {
    ContainerRef particles; 
    SettingsParam settings;
    try {
        cliParse(argc, argsv, particles, settings);
    } catch (YAMLReaderException& e) {
        SPDLOG_ERROR("YAML Reader failed with: {}", e.what());
        exit(-1);
    }
    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value(), settings.end_time.value());
    
    startSimulation(particles, settings);
    return 0;
}
