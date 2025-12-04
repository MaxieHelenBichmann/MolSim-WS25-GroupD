#include <spdlog/spdlog.h>

#include <cstdlib>
#include <memory>
#include <optional>
#include <stdexcept>

#include "exceptions/CLIException.h"
#include "exceptions/InputException.h"
#include "exceptions/MolSimException.h"
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

    // Phase 1: Parse CLI and read settings
    try {
        file_name = cliParse(argc, argsv);
        std::filesystem::path path = file_name;
        if (path.extension() == ".txt") {
            file_reader = std::make_unique<XVMReader>();
        } else if (path.extension() == ".yaml") {
            file_reader = std::make_unique<YAMLReader>();
        } else {
            SPDLOG_ERROR("Unsupported file extension: {}", path.extension().string());
            return EXIT_FAILURE;
        }
        file_reader->readSettings(settings, file_name);
        SPDLOG_INFO("Loaded settings from {}", file_name);
    } catch (const InputException&) {
        return EXIT_FAILURE;
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
            SPDLOG_DEBUG("Using gravitational force model");
            break;
        }
        case LENNARDJONES: {
            force = std::make_unique<LennardJonesForce>();
            SPDLOG_DEBUG("Using Lennard-Jones force model");
            break;
        }
    }

    // Phase 2: Read particles and run simulation
    try {
        if (settings.container_type == "SIMPLE") {
            SimpleContainer particle_container(settings.domain.getDimension(), settings.cutoff);
            file_reader->readParticles(particle_container, file_name);
            SPDLOG_INFO("Loaded {} particles from {}", particle_container.size(), file_name);
            SPDLOG_INFO("Simulation configured: {} particles, delta_t={}, t=[{}, {}]", particle_container.size(),
                        settings.delta_t, settings.start_time, settings.end_time);

            Simulation<SimpleContainer> simulation(particle_container, *force, settings, *writer);
            simulation.run();
        } else if (settings.container_type == "LINKED") {
            LinkedCellContainer particle_container{settings.domain.getDimension(), settings.cutoff};
            file_reader->readParticles(particle_container, file_name);
            SPDLOG_INFO("Loaded {} particles from {}", particle_container.size(), file_name);
            SPDLOG_INFO("Simulation configured: {} particles, delta_t={}, t=[{}, {}]", particle_container.size(),
                        settings.delta_t, settings.start_time, settings.end_time);

            Simulation<LinkedCellContainer> simulation(particle_container, *force, settings, *writer);
            simulation.run();
        } else {
            SPDLOG_ERROR("Unknown container type: {}", settings.container_type);
            return EXIT_FAILURE;
        }
    } catch (const InputException&) {
        SPDLOG_ERROR("Caught Input Exception, Exiting now!");
        return EXIT_FAILURE;
    } catch (const SimulationException&) {
        SPDLOG_ERROR("Caught Simulation Exception, Exiting now!");
        return EXIT_FAILURE;
    } catch (const MolSimException& e) {
        SPDLOG_ERROR("Unexpected MolSim error: {}, Exiting now!", e.what());
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
        SPDLOG_ERROR("Unexpected error: {}, Exiting now!", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
