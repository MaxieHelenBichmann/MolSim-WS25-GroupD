#include <benchmark/benchmark.h>

#include <memory>

#include "../code/simulationimpl/ForceCalcMethods.h"
#include "BenchmarkingUtils.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "physics/LennardJonesForce.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {
/**
 * @brief Benchmarks the unoptimized Force Calculation
 * Performs one calulateF() call on a container of 100-1000 particles
 *
 */
template <ForceConcept forceType>
void bmForceUnOptimized(benchmark::State& state) {
    SimpleContainer particles;
    size_t n = state.range(0);
    particles.reserve(n);

    for (size_t i = 0; i < n; i++) {
        particles.addParticle(randomParticle());
    }
    for ([[maybe_unused]] auto _ : state) {
        calculateFUnoptimized<SimpleContainer, forceType>(particles);
        benchmark::DoNotOptimize(particles);
    }
}

/**
 * @brief Benchmarks the optimized Force Calculation
 * Performs one calulateF() call on a container of 100-1000 particles
 *
 */
template <ForceConcept forceType>
void bmForceOptimized(benchmark::State& state) {
    SimpleContainer particles;
    size_t n = state.range(0);
    particles.reserve(n);

    for (size_t i = 0; i < n; i++) {
        particles.addParticle(randomParticle());
    }

    for ([[maybe_unused]] auto _ : state) {
        calculateF<SimpleContainer, forceType>(particles);
        benchmark::DoNotOptimize(particles);
    }
}

/**
 * @brief Benchmarks the alternative optimized Force Calculation
 * Performs one calulateF() call on a container of 100-1000 particles
 *
 */
template <ForceConcept forceType>
void bmForceOptimizedAlt(benchmark::State& state) {
    SimpleContainer particles;
    size_t n = state.range(0);
    particles.reserve(n);

    for (size_t i = 0; i < n; i++) {
        particles.addParticle(randomParticle());
    }

    for ([[maybe_unused]] auto _ : state) {
        calculateFAlt<SimpleContainer, forceType>(particles);
        benchmark::DoNotOptimize(particles);
    }
}

BENCHMARK(bmForceUnOptimized<LennardJonesForce>)
    ->RangeMultiplier(10)
    ->Range(100, 10000)
    ->Complexity()
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond)
    ->DisplayAggregatesOnly(true);
BENCHMARK(bmForceOptimized<LennardJonesForce>)
    ->RangeMultiplier(10)
    ->Range(100, 10000)
    ->Complexity()
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond)
    ->DisplayAggregatesOnly(true);
BENCHMARK(bmForceOptimizedAlt<LennardJonesForce>)
    ->RangeMultiplier(10)
    ->Range(100, 10000)
    ->Complexity()
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond)
    ->DisplayAggregatesOnly(true);
}  // namespace mol_sim
