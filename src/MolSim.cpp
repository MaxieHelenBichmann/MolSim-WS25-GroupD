#include <spdlog/spdlog.h>

#include "io/CLIParse.h"
#include "io/fileReader/YAMLReaderException.h"
#include "particles/container/SimpleContainer.h"
#include "physics/GravitationalForce.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

using namespace mol_sim;
// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char* argsv[]) {
    SimpleContainer particles;
    SettingsParam settings;
    try {
        cliParse(argc, argsv, particles, settings);
    } catch (YAMLReaderException& e) {
        SPDLOG_ERROR("YAML Reader failed with: {}", e.what());
        exit(-1);
    }
    settings.setDefaults();
    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value(), settings.end_time.value());  // NOLINT(bugprone-unchecked-optional-access)

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
