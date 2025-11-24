#include <benchmark/benchmark.h>

#include "particles/container/ContainerRef.h"
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
    LennardJonesForce lj_force;
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    size_t n = state.range(0);
    CuboidGenerator generator({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {n, n, n}, 1.0, 1.0, 0.5, 5.0, 1.0);
    SettingsParam settings(0.1, 0, 1000, 5.0, 1.0);
    Simulation<SimpleContainer, LennardJonesForce> simulation(part_container, lj_force, settings);
    for ([[maybe_unused]] auto _ : state) {
        generator.generateParticles(particles);
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
    Simulation<SimpleContainer, LennardJonesForce> simulation(part_container, lj_force, settings);
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
