/**
 * @file BoundaryBenchmark.cpp
 * @brief Benchmarks comparing different boundary condition implementations.
 *
 * Compares the performance of Reflecting (ghost particles) vs VelocityReflect
 * boundary conditions by measuring applyBoundary() calls in isolation.
 *
 * Filter: --benchmark_filter=Boundary/
 */
#include <benchmark/benchmark.h>

#include <array>
#include <memory>

#include "../code/boundaryimpl/PeriodicOld.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Periodic.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/boundaries/VelocityReflect.h"
#include "particles/container/LinkedCellContainer.h"
#include "physics/pairwiseforces/LennardJonesForce.h"
#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

/**
 * @brief Generates particles near boundaries for boundary benchmarking.
 * Places particles within sigma distance of each boundary face.
 */
void generateBoundaryParticles(LinkedCellContainer& container, R3 domain_size, size_t particles_per_face,
                               double sigma) {
    // Generate particles near each of the 6 boundary faces
    for (size_t face = 0; face < 6; ++face) {
        for (size_t i = 0; i < particles_per_face; ++i) {
            R3 pos;
            R3 vel = maxwellBoltzmannDistributedVelocity(1.0, 3);

            double margin = sigma * 0.5;
            double rand_x =
                margin + ((static_cast<double>(rand()) / RAND_MAX) * (domain_size[0] - 2 * margin));  // NOLINT
            double rand_y =
                margin + ((static_cast<double>(rand()) / RAND_MAX) * (domain_size[1] - 2 * margin));  // NOLINT
            double rand_z =
                margin + ((static_cast<double>(rand()) / RAND_MAX) * (domain_size[2] - 2 * margin));  // NOLINT

            switch (face) {
                case 0:  // LEFT (x = 0)
                    pos = {margin, rand_y, rand_z};
                    break;
                case 1:  // RIGHT (x = domain_size[0])
                    pos = {domain_size[0] - margin, rand_y, rand_z};
                    break;
                case 2:  // FRONT (y = 0)
                    pos = {rand_x, margin, rand_z};
                    break;
                case 3:  // BACK (y = domain_size[1])
                    pos = {rand_x, domain_size[1] - margin, rand_z};
                    break;
                case 4:  // LOWER (z = 0)
                    pos = {rand_x, rand_y, margin};
                    break;
                default:  // UPPER (z = domain_size[2])
                    pos = {rand_x, rand_y, domain_size[2] - margin};
                    break;
            }
            container.addParticle(pos, vel, 1.0, 5.0, sigma);
        }
    }
}

/**
 * @brief Benchmarks applyBoundary with ghost-particle Reflecting boundaries.
 * Measures time to apply boundary conditions to N boundary particles.
 */
void bmBoundaryReflecting(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    double sigma = 1.0;

    LinkedCellContainer container(domain_size, cutoff);
    generateBoundaryParticles(container, domain_size, n / 6, sigma);

    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Reflecting>(BoundaryLocation::LEFT, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::RIGHT, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::FRONT, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::BACK, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::UPPER, domain_size, false),
        std::make_unique<Reflecting>(BoundaryLocation::LOWER, domain_size, false)};

    LennardJonesForce force;

    for ([[maybe_unused]] auto _ : state) {
        for (auto it = container.boundaryBegin(); it != container.boundaryEnd(); ++it) {
            for (const auto& boundary : boundaries) {
                boundary->applyBoundary(*it, force);
            }
        }
        benchmark::ClobberMemory();
    }
}

/**
 * @brief Benchmarks applyBoundary with velocity-inversion VelocityReflect boundaries.
 * Measures time to apply boundary conditions to N boundary particles.
 */
void bmBoundaryVelocityReflect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    double sigma = 1.0;

    LinkedCellContainer container(domain_size, cutoff);
    generateBoundaryParticles(container, domain_size, n / 6, sigma);

    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<VelocityReflect>(BoundaryLocation::LEFT, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::RIGHT, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::FRONT, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::BACK, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::UPPER, domain_size),
        std::make_unique<VelocityReflect>(BoundaryLocation::LOWER, domain_size)};

    LennardJonesForce force;

    for ([[maybe_unused]] auto _ : state) {
        for (auto it = container.boundaryBegin(); it != container.boundaryEnd(); ++it) {
            for (const auto& boundary : boundaries) {
                boundary->applyBoundary(*it, force);
            }
        }
        benchmark::ClobberMemory();
    }
}
/**
 * @brief Benchmarks applyBoundary with velocity-inversion VelocityReflect boundaries.
 * Measures time to apply boundary conditions to N boundary particles.
 */
void bmBoundaryPeriodic(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    double sigma = 1.0;

    LinkedCellContainer container(domain_size, cutoff);
    generateBoundaryParticles(container, domain_size, n / 6, sigma);

    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<Periodic>(BoundaryLocation::LEFT, domain_size, cutoff, 3),
        std::make_unique<Periodic>(BoundaryLocation::RIGHT, domain_size, cutoff, 3),
        std::make_unique<Periodic>(BoundaryLocation::FRONT, domain_size, cutoff, 3),
        std::make_unique<Periodic>(BoundaryLocation::BACK, domain_size, cutoff, 3),
        std::make_unique<Periodic>(BoundaryLocation::UPPER, domain_size, cutoff, 3),
        std::make_unique<Periodic>(BoundaryLocation::LOWER, domain_size, cutoff, 3)};

    LennardJonesForce force;

    for ([[maybe_unused]] auto _ : state) {
        for (auto it = container.boundaryBegin(); it != container.boundaryEnd(); ++it) {
            for (const auto& boundary : boundaries) {
                boundary->applyBoundary(*it, force);
            }
        }
        benchmark::ClobberMemory();
    }
}

/**
 * @brief Benchmarks applyBoundary with velocity-inversion VelocityReflect boundaries.
 * Measures time to apply boundary conditions to N boundary particles.
 */
void bmBoundaryPeriodicOld(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    double sigma = 1.0;

    LinkedCellContainer container(domain_size, cutoff);
    generateBoundaryParticles(container, domain_size, n / 6, sigma);

    std::array<std::unique_ptr<Boundary>, 6> boundaries{
        std::make_unique<PeriodicOld>(BoundaryLocation::LEFT, domain_size, cutoff, 3),
        std::make_unique<PeriodicOld>(BoundaryLocation::RIGHT, domain_size, cutoff, 3),
        std::make_unique<PeriodicOld>(BoundaryLocation::FRONT, domain_size, cutoff, 3),
        std::make_unique<PeriodicOld>(BoundaryLocation::BACK, domain_size, cutoff, 3),
        std::make_unique<PeriodicOld>(BoundaryLocation::UPPER, domain_size, cutoff, 3),
        std::make_unique<PeriodicOld>(BoundaryLocation::LOWER, domain_size, cutoff, 3)};

    LennardJonesForce force;

    for ([[maybe_unused]] auto _ : state) {
        for (auto it = container.boundaryBegin(); it != container.boundaryEnd(); ++it) {
            for (const auto& boundary : boundaries) {
                boundary->applyBoundary(*it, force);
            }
        }
        benchmark::ClobberMemory();
    }
}

BENCHMARK(bmBoundaryReflecting)
    ->Name("Boundary/Reflecting")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmBoundaryVelocityReflect)
    ->Name("Boundary/VelocityReflect")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmBoundaryPeriodic)
    ->Name("Boundary/Periodic/Current")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);

BENCHMARK(bmBoundaryPeriodicOld)
    ->Name("Boundary/Periodic/Old")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
}  // namespace mol_sim
