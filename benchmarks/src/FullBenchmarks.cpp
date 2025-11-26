#include <benchmark/benchmark.h>

#include <memory>

#include "../code/linkedcellimpl/LinkedCellContainerDirect.h"
#include "../code/linkedcellimpl/LinkedCellContainerExplicit.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/container/SimpleContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {
/**
 * @brief Tests the simulation with the parameters given in Assignment 2
 * Particle counts are 40x8x1 + 8x8x1
 *
 */
static void bmSimulationGiven(benchmark::State& state) {
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator1({20.0, 0.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    SettingsParam settings(0.0005, 0, 20, 5.0, 1.0);
    settings.domain = {.dimension = {180.0, 90., 1.}};
    settings.setDefaults();
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<SimpleContainer> simulation(part_container, std::move(force_source), settings, std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(particles);
        generator2.generateParticles(particles);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}

void bmSimulationGivenCutOff(benchmark::State& state) {
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator1({20.0, 0.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    SettingsParam settings(0.0005, 0, 20, 5.0, 1.0);
    settings.cutoff = 3.0;
    settings.domain = {.dimension = {180.0, 90., 1.}};
    settings.setDefaults();
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<SimpleContainer> simulation(part_container, std::move(force_source), settings, std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(particles);
        generator2.generateParticles(particles);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}

void bmSimulationGivenLCDirect(benchmark::State& state) {
    LinkedCellContainerDirect part_container({180., 90., 1.}, 3.0);
    ContainerRef particles(part_container);
    CuboidGenerator generator1({20.0, 0.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    SettingsParam settings(0.0005, 0, 20, 5.0, 1.0);
    settings.cutoff = 3.0;
    settings.domain = {.dimension = {180.0, 90., 1.}};
    settings.setDefaults();
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<LinkedCellContainerDirect> simulation(part_container, std::move(force_source), settings,
                                                     std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(particles);
        generator2.generateParticles(particles);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}

void bmSimulationGivenLCExplicit(benchmark::State& state) {
    LinkedCellContainerExplicit part_container({180., 90., 1.}, 3.0);
    ContainerRef particles(part_container);
    CuboidGenerator generator1({20.0, 0.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    SettingsParam settings(0.0005, 0, 20, 5.0, 1.0);
    settings.cutoff = 3.0;
    settings.domain = {.dimension = {180.0, 90., 1.}};
    settings.setDefaults();
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<LinkedCellContainerExplicit> simulation(part_container, std::move(force_source), settings,
                                                       std::move(writer));
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(particles);
        generator2.generateParticles(particles);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}

BENCHMARK(bmSimulationGiven)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationGivenCutOff)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationGivenLCDirect)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);
BENCHMARK(bmSimulationGivenLCExplcit)->Unit(benchmark::kMillisecond)->Repetitions(5)->DisplayAggregatesOnly(true);

}  // namespace mol_sim
