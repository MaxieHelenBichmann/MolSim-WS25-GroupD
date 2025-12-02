#include <benchmark/benchmark.h>

#include <array>
#include <memory>

#include "../code/linkedcellimpl/LinkedCellContainerDirect.h"
#include "../code/linkedcellimpl/LinkedCellContainerExplicit.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/boundaries/Boundary.h"
#include "particles/container/SimpleContainer.h"
#include "physics/LennardJonesForce.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {

/**
 * @brief Generates a cuboid of particles directly into any ParticleContainer.
 */
template <ParticleContainer Container>
void generateCuboid(Container& particles, R3 position, R3 velocity, Vector<size_t, 3> num_particles, double mass,
                    double distance, double avg_velo, double epsilon, double sigma) {
    particles.reserve(num_particles[0] * num_particles[1] * num_particles[2]);
    for (size_t i = 0; i < num_particles[2]; i++) {
        for (size_t j = 0; j < num_particles[1]; j++) {
            for (size_t k = 0; k < num_particles[0]; k++) {
                R3 curr_pos = {position[0] + (static_cast<double>(k) * distance),
                               position[1] + (static_cast<double>(j) * distance),
                               position[2] + (static_cast<double>(i) * distance)};
                R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
                particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma);
            }
        }
    }
}
/**
 * @brief Tests the simulation with the parameters given in Assignment 2
 * Particle counts are 40x8x1 + 8x8x1
 *
 */
static void bmSimulationGiven(benchmark::State& state) {
    SimpleContainer part_container;
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 20;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;

    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT),  std::make_unique<Outflow>(BoundaryLocation::RIGHT),
        std::make_unique<Outflow>(BoundaryLocation::FRONT), std::make_unique<Outflow>(BoundaryLocation::BACK),
        std::make_unique<Outflow>(BoundaryLocation::UPPER), std::make_unique<Outflow>(BoundaryLocation::LOWER)};
    settings.domain = Domain({180.0, 90., 1.}, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<SimpleContainer> simulation(part_container, std::move(force_source), settings, std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

void bmSimulationGivenCutOff(benchmark::State& state) {
    SimpleContainer part_container;
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 20;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT),  std::make_unique<Outflow>(BoundaryLocation::RIGHT),
        std::make_unique<Outflow>(BoundaryLocation::FRONT), std::make_unique<Outflow>(BoundaryLocation::BACK),
        std::make_unique<Outflow>(BoundaryLocation::UPPER), std::make_unique<Outflow>(BoundaryLocation::LOWER)};
    settings.domain = Domain({180.0, 90., 1.}, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<SimpleContainer> simulation(part_container, std::move(force_source), settings, std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

void bmSimulationGivenLCDirect(benchmark::State& state) {
    LinkedCellContainerDirect part_container({180., 90., 1.}, 3.0);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 20;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT),  std::make_unique<Outflow>(BoundaryLocation::RIGHT),
        std::make_unique<Outflow>(BoundaryLocation::FRONT), std::make_unique<Outflow>(BoundaryLocation::BACK),
        std::make_unique<Outflow>(BoundaryLocation::UPPER), std::make_unique<Outflow>(BoundaryLocation::LOWER)};
    settings.domain = Domain({180.0, 90., 1.}, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<LinkedCellContainerDirect> simulation(part_container, std::move(force_source), settings,
                                                     std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

void bmSimulationGivenLCExplicit(benchmark::State& state) {
    LinkedCellContainerExplicit part_container({180., 90., 1.}, 3.0);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 20;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;

    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT),  std::make_unique<Outflow>(BoundaryLocation::RIGHT),
        std::make_unique<Outflow>(BoundaryLocation::FRONT), std::make_unique<Outflow>(BoundaryLocation::BACK),
        std::make_unique<Outflow>(BoundaryLocation::UPPER), std::make_unique<Outflow>(BoundaryLocation::LOWER)};
    settings.domain = Domain({180.0, 90., 1.}, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<LinkedCellContainerExplicit> simulation(part_container, std::move(force_source), settings,
                                                       std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

BENCHMARK(bmSimulationGiven)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationGivenCutOff)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationGivenLCDirect)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);
BENCHMARK(bmSimulationGivenLCExplicit)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);

}  // namespace mol_sim
