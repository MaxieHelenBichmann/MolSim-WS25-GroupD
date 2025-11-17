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
    for ([[maybe_unused]] auto _ : state) {
        SimpleContainer part_container;
        CuboidGenerator generator1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {40U, 8U, 1U}, 1.0, 1.0, 0.1, 5.0, 1.0);
        CuboidGenerator generator2({15.0, 15.0, 0.0}, {0.0, -10.0, 0.0}, {8U, 8U, 1U}, 1.0, 1.0, 0.1, 5.0, 1.0);
        generator1.generateParticles(part_container);
        generator2.generateParticles(part_container);
        SettingsParam settings(0.14, 0, 1000, 0.0, 1.0);
        auto force_source = std::make_unique<LennardJonesForce>();
        auto writer = std::make_unique<XYZWriter>();
        Simulation<SimpleContainer> simulation(part_container, std::move(force_source), settings, std::move(writer));
        simulation.run();
    }
}
BENCHMARK(bmSimulationGiven)->Repetitions(10)->Unit(benchmark::kMillisecond);
}  // namespace mol_sim
