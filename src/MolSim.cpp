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
    SimpleContainer particles; /** @deprecated This should become obsolete. Still need to adjust XVM reader tho. YAMLReader already done. */
    SettingsParam settings;
    try {
        cliParse(argc, argsv, particles, settings);
    } catch (YAMLReaderException& e) {
        SPDLOG_ERROR("YAML Reader failed with: {}", e.what());
        exit(-1);
    }
    SPDLOG_INFO("Simulation configured with {} particles, delta_t={} end_time={}", particles.size(),
                settings.delta_t.value(), settings.end_time.value());

    /** TODO: clean this up
     * no idea how tho. if we were using normal super classes instead of the 2 concepts this wouldn't be an issue.
     * the thing is that we determine both of the parameters at runtime meaning theres nothing we can really do in terms of
     * compile time programming. also idk how much std::variant (in combination with std::visit) would really help here, just
     * in terms of maintainability (and code aesthetics ig (tho that parts less important)).
     */
    if (settings.domain_type.value().compare("SIMPLE") == 0) {
    switch (settings.force.value()) {
        case GRAVITATIONAL: {
            GravitationalForce grav_force;
            Simulation<SimpleContainer, GravitationalForce> simulation(settings, grav_force);
            simulation.run();
            return 0;
        }
        case LENNARDJONES: {
            LennardJonesForce lj_force;
            Simulation<SimpleContainer, LennardJonesForce> simulation(settings, lj_force);
            simulation.run();
            return 0;
        }
        default:
            return 0;
    }
    } else if (settings.domain_type.value().compare("LINKED") == 0) {
    switch (settings.force.value()) {
        case GRAVITATIONAL: {
            GravitationalForce grav_force;
            Simulation<LinkedCellContainer, GravitationalForce> simulation(settings, grav_force);
            simulation.run();
            return 0;
        }
        case LENNARDJONES: {
            LennardJonesForce lj_force;
            Simulation<LinkedCellContainer, LennardJonesForce> simulation(settings, lj_force);
            simulation.run();
            return 0;
        }
        default:
            return 0;

    }
}
