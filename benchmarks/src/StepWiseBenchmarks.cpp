/**
 * @file StepWiseBenchmarks.cpp
 * @brief Benchmarks for single simulation step performance comparison.
 *
 * Compares LinkedCellContainer vs SimpleContainer for a single simulation step,
 * useful for measuring per-step overhead without full simulation time.
 *
 * Filter: --benchmark_filter=Simulation/SingleStep/
 */
#include <benchmark/benchmark.h>

#include <memory>

#include "io/checkpointWriter/XVMWriterCP.h"
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
 * @brief Benchmarks single simulation step with LinkedCellContainer.
 * Runs one step with NxNxN cuboid (N = 10-20) particles.
 */
void bmSimulationSingleStepLinkedCell(benchmark::State& state) {
    size_t n = state.range(0);
    LinkedCellContainer part_container({180., 90., 1.}, 3.0);
    ContainerRef particles(part_container);

    CuboidGenerator generator({60.0, 60.0, 0.0}, {0.0, 0.0, 0.0}, {n, n, n}, 1.0, 1.1225, 0.1, 5.0, 1.0, 0.01);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 0.0005;
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
    auto cp_writer = std::make_unique<XVMWriterCP>();
    auto stat_writer = std::make_unique<StatsWriter>();
    generator.generateParticles(particles);

    Simulation<LinkedCellContainer> simulation(part_container, *force_source, settings, *writer, *cp_writer,
                                               *stat_writer);
    for ([[maybe_unused]] auto _ : state) {
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}

/**
 * @brief Benchmarks single simulation step with SimpleContainer (direct sum).
 * Runs one step with NxNxN cuboid (N = 10-20) particles.
 */
void bmSimulationSingleStepDirectSum(benchmark::State& state) {
    size_t n = state.range(0);
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator({60.0, 60.0, 0.0}, {0., 0., 0.}, {n, n, n}, 1.0, 1.1225, 0.1, 5.0, 1.0, 0.01);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 0.0005;
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
    auto cp_writer = std::make_unique<XVMWriterCP>();
    auto stat_writer = std::make_unique<StatsWriter>();
    generator.generateParticles(particles);
    Simulation<SimpleContainer> simulation(part_container, *force_source, settings, *writer, *cp_writer, *stat_writer);
    for ([[maybe_unused]] auto _ : state) {
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(particles);
    }
}
BENCHMARK(bmSimulationSingleStepLinkedCell)
    ->Name("Simulation/SingleStep/LinkedCell")
    ->DenseRange(10, 20, 2)
    ->Unit(benchmark::kMillisecond)
    ->Complexity()
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationSingleStepDirectSum)
    ->Name("Simulation/SingleStep/DirectSum")
    ->DenseRange(10, 20, 2)
    ->Unit(benchmark::kMillisecond)
    ->Complexity()
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
}  // namespace mol_sim
