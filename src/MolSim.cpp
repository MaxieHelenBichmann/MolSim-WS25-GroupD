#include <spdlog/spdlog.h>

#include "io/CLIParse.h"
#include "io/fileReader/YAMLReader.h"
#include "io/fileReader/YAMLReaderException.h"
#include "particles/container/SimpleContainer.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

using namespace mol_sim;
// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char* argsv[]) {
    YAMLReader file_reader;
    SimpleContainer particles;
    SettingsParam settings;
    try {
        cliParse(argc, argsv, file_reader, particles, settings);
    } catch (YAMLReaderException& e) {
        SPDLOG_ERROR("YAML Reader failed with: {}", e.what());
        exit(-1);
    }
    settings.setDefaults();

    // explicit .value_or despite call to setDefaults to make warnings disappear
    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value(), settings.end_time.value());  // NOLINT((bugprone-unchecked-optional-access))
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.run();
}
