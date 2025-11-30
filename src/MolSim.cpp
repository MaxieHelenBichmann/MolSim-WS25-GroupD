#include <spdlog/spdlog.h>

#include <stdexcept>

#include "io/CLIParse.h"
#include "io/FileReader.h"
#include "io/fileReader/YAMLReader.h"
#include "io/fileReader/YAMLReaderException.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/SimpleContainer.h"
#include "physics/GravitationalForce.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

using namespace mol_sim;

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char* argsv[]) {
    SettingsParam settings;
    std::unique_ptr<FileReader> file_reader;
    std::string file_name;
    try {
        file_name = cliParse(argc, argsv, settings);
        std::filesystem::path path = file_name;
        if (path.extension() == ".txt") {
            file_reader = std::make_unique<XVMReader>();
        } else if (path.extension() == ".yaml") {
            file_reader = std::make_unique<YAMLReader>();
        } else {
            SPDLOG_ERROR("Unsupported File Extension!");
            exit(-1);
        }
        file_reader->readSettings(settings, file_name);
        settings.setDefaults();
    } catch (std::runtime_error& e) {
        SPDLOG_ERROR("Settings parser failed with: {}", e.what());
        exit(-1);
    }
    if (settings.container_type.value() == "SIMPLE") {
        SimpleContainer particle_container(settings.domain.value().getDimension(), settings.cutoff.value());
        try {
            file_reader->readParticles(particle_container, file_name);
        } catch (std::runtime_error& e) {
            SPDLOG_ERROR("Reader failed with: {}", e.what());
            exit(-1);
        }
        SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particle_container.size(),
                    settings.delta_t.value(), settings.end_time.value());
        LennardJonesForce f;
        Simulation<SimpleContainer, LennardJonesForce> simulation(particle_container, f, settings);
        simulation.run();
    } else if (settings.container_type.value() == "LINKED") {
        LinkedCellContainer particle_container{settings.domain.value().getDimension(), settings.cutoff.value()};
        try {
            file_reader->readParticles(particle_container, file_name);
        } catch (std::runtime_error& e) {
            SPDLOG_ERROR("Reader failed with: {}", e.what());
            exit(-1);
        }
        SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particle_container.size(),
                    settings.delta_t.value(), settings.end_time.value());
        LennardJonesForce f;
        Simulation<LinkedCellContainer, LennardJonesForce> simulation(particle_container, f, settings);
        simulation.run();
    }
}
