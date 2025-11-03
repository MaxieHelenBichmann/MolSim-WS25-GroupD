#include <spdlog/spdlog.h>

#include "io/CLIParse.h"
#include "io/fileReader/YAMLReader.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

using namespace mol_sim;

int main(int argc, char* argsv[]) {
    YAMLReader file_reader;
    SimpleContainer particles;
    SettingsParam settings;
    cliParse(argc, argsv, file_reader, particles, settings);
    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(), settings.delta_t,
                settings.end_time);
    Simulation<SimpleContainer> simulation(particles, GRAVITATIONAL, settings);
    simulation.run();
}
