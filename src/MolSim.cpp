#include <spdlog/spdlog.h>

#include "io/CLIParse.h"
#include "io/fileReader/YAMLReader.h"
#include "particles/container/SimpleContainer.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"
#include "utils/Defaults.h"

using namespace mol_sim;

int main(int argc, char* argsv[]) {
    YAMLReader file_reader;
    SimpleContainer particles;
    SettingsParam settings;
    cliParse(argc, argsv, file_reader, particles, settings);
    setDefaults(settings);

    //explicit .value_or despite call to setDefaults to make warnings disappear
    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value_or(DELTA_T_DEFAULT), settings.end_time.value_or(END_TIME_DEFAULT));
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.run();
}
