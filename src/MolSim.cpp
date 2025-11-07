#include <spdlog/spdlog.h>

#include "io/CLIParse.h"
#include "io/fileReader/YAMLReaderException.h"
#include "particles/container/SimpleContainer.h"
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
                settings.delta_t.value(), settings.end_time.value());  // NOLINT((bugprone-unchecked-optional-access))
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, settings);
    simulation.run();
}
