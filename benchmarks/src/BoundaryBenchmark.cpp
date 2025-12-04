/**
 * @file BoundaryBenchmark.cpp
 * @brief Benchmarks comparing different boundary condition implementations.
 *
 * Compares the performance of Reflecting (ghost particles) vs VelocityReflect
 * boundary conditions in a full simulation with Lennard-Jones interactions.
 *
 * Filter: --benchmark_filter=Boundary/
 */
#include <benchmark/benchmark.h>

#include <array>
#include <memory>

#include "io/outputWriter/XYZWriter.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/boundaries/VelocityReflect.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {

/**
 * @brief Benchmarks simulation with ghost-particle Reflecting boundaries.
 * Runs a full 20s simulation with 2400 particles (100x20x1 + 20x20x1 cuboids).
 */
void bmBoundaryReflecting(benchmark::State& state) {
    LinkedCellContainer part_container({180., 90., 1.}, 3.0);
    ContainerRef ref(part_container);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 5.;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;

    R3 domain_size = {180.0, 90., 1.};
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Reflecting>(BoundaryLocation::LEFT, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::RIGHT, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::FRONT, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::BACK, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::UPPER, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::LOWER, domain_size, false)};
    settings.domain = Domain(domain_size, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    CuboidGenerator generator1({20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    Simulation<LinkedCellContainer> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(ref);
        generator2.generateParticles(ref);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

/**
 * @brief Benchmarks simulation with velocity-inversion Reflecting boundaries.
 * Runs a full 20s simulation with 2400 particles (100x20x1 + 20x20x1 cuboids).
 */
void bmBoundaryVelocityReflect(benchmark::State& state) {
    LinkedCellContainer part_container({180., 90., 1.}, 3.0);
    ContainerRef ref(part_container);
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0;
    settings.end_time = 5.;
    settings.epsilon = 5.0;
    settings.sigma = 1.0;
    settings.cutoff = 3.0;

    R3 domain_size = {180.0, 90., 1.};
    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<VelocityReflect>(BoundaryLocation::LEFT, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::RIGHT, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::FRONT, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::BACK, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::UPPER, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::LOWER, domain_size)};
    settings.domain = Domain(domain_size, std::move(boundaries));
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    CuboidGenerator generator1({20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    CuboidGenerator generator2({70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
    Simulation<LinkedCellContainer> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        generator1.generateParticles(ref);
        generator2.generateParticles(ref);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

BENCHMARK(bmBoundaryReflecting)
    ->Name("Boundary/Reflecting")
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK(bmBoundaryVelocityReflect)
    ->Name("Boundary/VelocityReflect")
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

}  // namespace mol_sim
