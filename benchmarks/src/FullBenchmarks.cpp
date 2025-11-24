#include <benchmark/benchmark.h>

#include <memory>

#include "io/outputWriter/XYZWriter.h"
#include "particles/container/SimpleContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"
namespace mol_sim {
/**
 * @brief Tests the simulation function with a large number of particles.
 * Particles range from 2x2x2 to 64x64x64
 *
 */
[[maybe_unused]] static void bmSimulationBig(benchmark::State& state) {
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        SimpleContainer part_container;
        CuboidGenerator generator({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {n, n, n}, 1.0, 1.0, 0.5, 5.0, 1.0);
        generator.generateParticles(part_container);
        SettingsParam settings(0.1, 0, 1000, 5.0, 1.0);
        auto force_source = std::make_unique<LennardJonesForce>();
        auto writer = std::make_unique<XYZWriter>();
        Simulation<SimpleContainer> simulation(part_container, std::move(force_source), settings, std::move(writer));
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}
BENCHMARK(bmSimulationBig)->RangeMultiplier(2)->Range(2, 2 << 6)->Repetitions(10)->Unit(benchmark::kMillisecond);
/**
 * @brief Tests the simulation with the parameters given in Assignment 2
 * Particle counts are 40x8x1 + 8x8x1
 *
 */
static void bmSimulationGiven(benchmark::State& state) {
    LennardJonesForce lj_force;
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator1({20.0, 0.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {10U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    SettingsParam settings(0.0005, 0, 20, 5.0, 1.0);
    settings.cutoff = 3.0;
    settings.domain = {.dimension = {180.0, 9.0, 1.}};
    settings.setDefaults();
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<SimpleContainer> simulation(part_container, std::move(force_source), settings, std::move(writer));
    simulation.run();
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(particles);
        generator2.generateParticles(particles);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}

void bmSimulationGivenCutOff(benchmark::State& state) {}

void bmSimulationGivenLC(benchmark::State& state) {}

BENCHMARK(bmSimulationGiven)->Repetitions(10)->Unit(benchmark::kMillisecond);
}  // namespace mol_sim
