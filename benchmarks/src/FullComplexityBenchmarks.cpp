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

#include <cmath>
#include <cstddef>
#include <memory>

#include "io/checkpointWriter/XVMWriterCP.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/boundaries/Outflow.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/SimpleContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/pairwiseforces/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {

/**
 * @brief Benchmarks LinkedCellContainer simulation complexity.
 * Runs simulation with NxNx1 cuboid to measure O(n) scaling.
 * Parameter is particle count (must be perfect square).
 */
void bmSimulationComplexityLinkedCell(benchmark::State& state) {
    size_t particle_count = state.range(0);
    auto n = static_cast<size_t>(std::sqrt(particle_count));
    double spacing = 1.1225;
    double grid_size = static_cast<double>(n) * spacing;
    double margin = 100.0;
    R3 domain_size = {grid_size + (2 * margin), grid_size + (2 * margin), 1.};
    LinkedCellContainer part_container(domain_size, 3.0);
    ContainerRef particles(part_container);

    CuboidGenerator generator({margin, margin, 0.0}, {0.0, 0.0, 0.0}, {n, n, static_cast<size_t>(1)}, 1.0, spacing, 0.3,
                              5.0, 1.0, 0.01);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 1.0;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::RIGHT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::FRONT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::BACK, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::UPPER, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::LOWER, domain_size)};
    settings.domain = Domain(domain_size, std::move(boundaries));
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    auto cp_writer = std::make_unique<XVMWriterCP>();
    auto writer = std::make_unique<XYZWriter>();

    Simulation<LinkedCellContainer> simulation(part_container, pairwise_forces, single_forces, settings, *writer,
                                               *cp_writer);
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
 * Runs simulation with NxNx1 cuboid to measure O(n²) scaling.
 * Parameter is particle count (must be perfect square).
 */
void bmSimulationComplexityDirectSum(benchmark::State& state) {
    size_t particle_count = state.range(0);
    auto n = static_cast<size_t>(std::sqrt(particle_count));
    double spacing = 1.1225;
    double grid_size = static_cast<double>(n) * spacing;
    double margin = 10.0;
    R3 domain_size = {grid_size + (2 * margin), grid_size + (2 * margin), 1.};
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator({margin, margin, 0.0}, {0., 0., 0.}, {n, n, static_cast<size_t>(1)}, 1.0, spacing, 0.3,
                              5.0, 1.0, 0.01);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 1.0;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Outflow>(BoundaryLocation::LEFT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::RIGHT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::FRONT, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::BACK, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::UPPER, domain_size),
        std::make_unique<Outflow>(BoundaryLocation::LOWER, domain_size)};
    settings.domain = Domain(domain_size, std::move(boundaries));
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    auto writer = std::make_unique<XYZWriter>();
    auto cp_writer = std::make_unique<XVMWriterCP>();
    Simulation<SimpleContainer> simulation(part_container, pairwise_forces, single_forces, settings, *writer,
                                           *cp_writer);
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
    ->Arg(100)
    ->Arg(400)
    ->Arg(900)
    ->Arg(1600)
    ->Arg(2500)
    //->Arg(3600)
    //->Arg(4900)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oN)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationComplexityDirectSum)
    ->Name("Simulation/Complexity/DirectSum")
    ->Arg(100)
    ->Arg(400)
    ->Arg(900)
    ->Arg(1600)
    ->Arg(2500)
    //->Arg(3600)
    //->Arg(4900)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oNSquared)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
}  // namespace mol_sim
