#include <spdlog/spdlog.h>

#include <memory>
#include <optional>

#include "exceptions/CLIException.h"
#include "exceptions/SimulationException.h"
#include "io/CLIParse.h"
#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
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
    } catch (const CLIException& e) {
        SPDLOG_ERROR("CLI parsing failed: {}", e.what());
        exit(-1);
    }

    settings.setDefaults();

    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value(), settings.end_time.value());

    std::unique_ptr<OutputWriter> writer;
#ifdef ENABLE_VTK_OUTPUT
    writer = std::make_unique<VTKWriter>();
#else
    writer = std::make_unique<XYZWriter>();
#endif

    std::unique_ptr<ForceSource> force;
    switch (settings.force.value()) {
        case GRAVITATIONAL: {
            force = std::make_unique<GravitationalForce>();
            break;
        }
        case LENNARDJONES: {
            force = std::make_unique<LennardJonesForce>();
            break;
        }
    }
    try {
        Simulation<SimpleContainer> simulation(particles, std::move(force), settings, std::move(writer));
        simulation.run();
    } catch (SimulationException& e) {
        exit(-1);
    }

    return 0;
}
