#include <spdlog/spdlog.h>

#include "io/CLIParse.h"
#include "io/fileReader/YAMLReader.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

using namespace mol_sim;
/**
 * @brief Default value for the delta_t parameter of the simulation.
 *
 */
constexpr double DELTA_T_DEFAULT = 0.014;
/**
 * @brief Default value for the end_time parameter of the simulation.
 *
 */
constexpr double END_TIME_DEFAULT = 1000;
/**
 * @brief Default value for the start_time parameter of the simulation.
 *
 */
constexpr double START_TIME_DEFAULT = 0;

int main(int argc, char* argsv[]) {
    YAMLReader file_reader;
    SimpleContainer particles;
    SettingsParam settings;
    cliParse(argc, argsv, file_reader, particles, settings);

    // Provide default values for settings that have not been set
    if (!settings.delta_t.has_value()) {
        settings.delta_t = DELTA_T_DEFAULT;
    }
    if (!settings.end_time.has_value()) {
        settings.end_time = END_TIME_DEFAULT;
    }
    if (!settings.start_time.has_value()) {
        settings.start_time = START_TIME_DEFAULT;
    }

    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value(), settings.end_time.value());
    Simulation<SimpleContainer> simulation(particles, GRAVITATIONAL, settings);
    simulation.run();
}
