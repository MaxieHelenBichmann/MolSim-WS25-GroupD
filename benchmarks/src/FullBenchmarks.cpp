/**
 * @file FullBenchmarks.cpp
 * @brief End-to-end simulation benchmarks comparing different container implementations.
 *
 * Runs full simulations with Assignment 3 Task 3 parameters to compare performance
 * of SimpleContainer, LinkedCellContainerDirect, and LinkedCellContainerExplicit.
 *
 * Filter: --benchmark_filter=Simulation/Full/
 */
#include <benchmark/benchmark.h>

#include <array>
#include <memory>

#include "../code/linkedcellimpl/LinkedCellContainerDirect.h"
#include "../code/linkedcellimpl/LinkedCellContainerExplicit.h"
#include "../code/linkedcellimpl/SimpleContainerControl.h"
#include "../code/simulationimpl/SimulationBenchmark.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Outflow.h"
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
 * @brief Benchmarks full simulation with SimpleContainer using Assignment 3 parameters.
 * Particle counts are 100x20x1 + 20x20x1 (2400 total), 20s simulation time.
 */
static void bmSimulationFullSimple(benchmark::State& state) {
    SimpleContainerControl part_container;
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
    SimulationBenchmark<SimpleContainerControl> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        part_container.clear();
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

/**
 * @brief Benchmarks full simulation with SimpleContainer using cutoff radius.
 * Particle counts are 100x20x1 + 20x20x1 (2400 total), 20s simulation time.
 */
void bmSimulationFullSimpleCutoff(benchmark::State& state) {
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
    SimpleContainerControl part_container(domain_size, settings.cutoff);
    SimulationBenchmark<SimpleContainerControl> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        part_container.clear();
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

/**
 * @brief Benchmarks full simulation with LinkedCellContainerDirect.
 * Particle counts are 100x20x1 + 20x20x1 (2400 total), 1s simulation time.
 */
void bmSimulationFullLinkedCellDirect(benchmark::State& state) {
    LinkedCellContainerDirect part_container({180., 90., 1.}, 3.0);
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
    SimulationBenchmark<LinkedCellContainerDirect> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        part_container.clear();
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

/**
 * @brief Benchmarks full simulation with LinkedCellContainerExplicit.
 * Particle counts are 100x20x1 + 20x20x1 (2400 total), 1s simulation time.
 */
void bmSimulationFullLinkedCellExplicit(benchmark::State& state) {
    LinkedCellContainerExplicit part_container({180., 90., 1.}, 3.0);
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
    SimulationBenchmark<LinkedCellContainerExplicit> simulation(part_container, *force_source, settings, *writer);
    for ([[maybe_unused]] auto _ : state) {
        part_container.clear();
        generateCuboid(part_container, {20.0, 20.0, 0.0}, {0., 0.0, 0.0}, {100U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0, 1.0);
        generateCuboid(part_container, {70.0, 60.0, 0.0}, {0.0, -10.0, 0.0}, {20U, 20U, 1U}, 1.0, 1.1225, 0.1, 5.0,
                       1.0);
        benchmark::ClobberMemory();
        simulation.run();
        benchmark::DoNotOptimize(part_container);
    }
}

BENCHMARK(bmSimulationFullSimple)
    ->Name("Simulation/Full/Simple")
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationFullSimpleCutoff)
    ->Name("Simulation/Full/SimpleCutoff")
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmSimulationFullLinkedCellDirect)
    ->Name("Simulation/Full/LinkedCellDirect")
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK(bmSimulationFullLinkedCellExplicit)
    ->Name("Simulation/Full/LinkedCellExplicit")
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

}  // namespace mol_sim
