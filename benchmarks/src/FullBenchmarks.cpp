#include <benchmark/benchmark.h>

#include "BenchmarkingUtils.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/ForceSource.h"
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
    setUpLogging();
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    size_t n = state.range(0);
    CuboidGenerator generator({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {n, n, n}, 1.0, 1.0, 0.5, 5.0, 1.0);
    SettingsParam settings(0.1, 0, 1000, LENNARDJONES, 5.0, 1.0);
    Simulation<SimpleContainer, LennardJonesForce> simulation(part_container, settings);
    for ([[maybe_unused]] auto _ : state) {
        generator.generateParticles(particles);
        simulation.run();
    }
}
BENCHMARK(bmSimulationBig)->RangeMultiplier(2)->Range(2, 2 << 6)->Repetitions(10)->Unit(benchmark::kMillisecond);
/**
 * @brief Tests the simulation with the parameters given in Assignment 2
 * Particle counts are 40x8x1 + 8x8x1
 *
 */
static void bmSimulationGiven(benchmark::State& state) {
    setUpLogging();
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {40U, 8U, 1U}, 1.0, 1.0, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({15.0, 15.0, 0.0}, {0.0, -10.0, 0.0}, {8U, 8U, 1U}, 1.0, 1.0, 0.1, 5.0, 1.0);
    SettingsParam settings(0.14, 0, 1000, LENNARDJONES, 0.0, 1.0);
    Simulation<SimpleContainer, LennardJonesForce> simulation(part_container, settings);
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(particles);
        generator2.generateParticles(particles);
        simulation.run();
    }
}
BENCHMARK(bmSimulationGiven)->Repetitions(10)->Unit(benchmark::kMillisecond);
}  // namespace mol_sim
