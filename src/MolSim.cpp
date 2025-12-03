#include <spdlog/spdlog.h>

#include <memory>
#include <optional>
#include <stdexcept>

#include "exceptions/CLIException.h"
#include "exceptions/SimulationException.h"
#include "exceptions/YAMLReaderException.h"
#include "io/CLIParse.h"
#include "io/FileReader.h"
#include "io/fileReader/YAMLReader.h"
#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
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
    // Phase 1 reading in settings
    try {
        file_name = cliParse(argc, argsv);
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
    } catch (const CLIException& e) {
        SPDLOG_ERROR("Settings parser failed with: {}", e.what());
        exit(-1);
    }

    std::unique_ptr<OutputWriter> writer;
#ifdef ENABLE_VTK_OUTPUT
    writer = std::make_unique<VTKWriter>();
#else
    writer = std::make_unique<XYZWriter>();
#endif

    std::unique_ptr<ForceSource> force;

    switch (settings.force) {
        case GRAVITATIONAL: {
            force = std::make_unique<GravitationalForce>();
            break;
        }
        case LENNARDJONES: {
            force = std::make_unique<LennardJonesForce>();
            break;
        }
    }

    // Phase 2 reading in Particles + running Simulation
    if (settings.container_type == "SIMPLE") {
        SimpleContainer particle_container(settings.domain.getDimension(), settings.cutoff);
        try {
            file_reader->readParticles(particle_container, file_name);
        } catch (std::runtime_error& e) {
            SPDLOG_ERROR("Reader failed with: {}", e.what());
            exit(-1);
        }
        SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particle_container.size(),
                    settings.delta_t, settings.end_time);

        try {
            Simulation<SimpleContainer> simulation(particle_container, *force, settings, *writer);
            simulation.run();
        } catch (SimulationException& e) {
            exit(-1);
        };
    } else if (settings.container_type == "LINKED") {
        LinkedCellContainer particle_container{settings.domain.getDimension(), settings.cutoff};
        try {
            file_reader->readParticles(particle_container, file_name);
        } catch (std::runtime_error& e) {
            SPDLOG_ERROR("Reader failed with: {}", e.what());
            exit(-1);
        }
        SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particle_container.size(),
                    settings.delta_t, settings.end_time);
        try {
            Simulation<LinkedCellContainer> simulation(particle_container, *force, settings, *writer);
            simulation.run();
        } catch (SimulationException& e) {
            exit(-1);
        }
    }
}
