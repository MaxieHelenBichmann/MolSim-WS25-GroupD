#include <spdlog/spdlog.h>

#include "exceptions/CLIException.h"
#include "io/CLIParse.h"
#include "particles/container/SimpleContainer.h"
#include "physics/GravitationalForce.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

using namespace mol_sim;

int main(int argc, char* argsv[]) {
    SimpleContainer particles;
    SettingsParam settings;
    try {
        cliParse(argc, argsv, particles, settings);
    } catch (const CLIException& e) {
        SPDLOG_ERROR("CLI parsing failed: {}", e.what());
        exit(-1);
    }
    settings.setDefaults();
    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value(), settings.end_time.value());

    switch (settings.force.value()) {
        case GRAVITATIONAL: {
            GravitationalForce grav_force;
            Simulation<SimpleContainer, GravitationalForce> simulation(particles, grav_force, settings);
            simulation.run();
            return 0;
        }
        case LENNARDJONES: {
            LennardJonesForce lj_force;
            Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
            simulation.run();
            return 0;
        }
        default:
            return 0;
    }
}
