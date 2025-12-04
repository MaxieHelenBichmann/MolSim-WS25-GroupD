/**
 * @file FullComplexityBenchmarks.cpp
 * @brief Benchmarks for measuring simulation complexity scaling (O(n) vs O(n²)).
 *
 * Compares LinkedCellContainer (expected O(n)) vs SimpleContainer (expected O(n²))
 * across varying particle counts to validate algorithmic complexity.
 *
 * Filter: --benchmark_filter=Simulation/Complexity/
 */
#include <benchmark/benchmark.h>

#include <cstddef>
#include <memory>

#include "io/outputWriter/XYZWriter.h"
#include "particles/boundaries/Outflow.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/SimpleContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {

/**
 * @brief Benchmarks LinkedCellContainer simulation complexity.
 * Runs 1s simulation with NxNxN cuboid (N = 10-20) to measure O(n) scaling.
 */
void bmSimulationComplexityLinkedCell(benchmark::State& state) {
    size_t n = state.range(0);
    LinkedCellContainer part_container({180., 90., 1.}, 3.0);
    ContainerRef particles(part_container);

    CuboidGenerator generator({60.0, 60.0, 0.0}, {0.0, 0.0, 0.0}, {n, n, static_cast<size_t>(1)}, 1.0, 1.1225, 0.3, 5.0,
                              1.0);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 1.0;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;
    R3 domain_size = {180.0, 90., 1.};
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::RIGHT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::FRONT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::BACK, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::UPPER, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::LOWER, domain_size)};
    settings.domain = Domain(domain_size, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();

    Simulation<LinkedCellContainer> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        particles.clear();
        generator.generateParticles(particles);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}

/**
 * @brief Benchmarks SimpleContainer simulation complexity.
 * Runs 2s simulation with NxNxN cuboid (N = 10-20) to measure O(n²) scaling.
 */
void bmSimulationComplexityDirectSum(benchmark::State& state) {
    size_t n = state.range(0);
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator({60.0, 60.0, 0.0}, {0., 0., 0.}, {n, n, static_cast<size_t>(1)}, 1.0, 1.1225, 0.3, 5.0,
                              1.0);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 1.0;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    R3 domain_size = {180.0, 90., 1.};
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::RIGHT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::FRONT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::BACK, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::UPPER, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::LOWER, domain_size)};
    settings.domain = Domain(domain_size, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<SimpleContainer> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        particles.clear();
        generator.generateParticles(particles);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}
BENCHMARK(bmSimulationComplexityLinkedCell)
    ->Name("Simulation/Complexity/LinkedCell")
    ->DenseRange(10, 20, 2)
    ->Unit(benchmark::kMillisecond)
    ->Complexity()
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationComplexityDirectSum)
    ->Name("Simulation/Complexity/DirectSum")
    ->DenseRange(10, 20, 2)
    ->Unit(benchmark::kMillisecond)
    ->Complexity()
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
}  // namespace mol_sim
