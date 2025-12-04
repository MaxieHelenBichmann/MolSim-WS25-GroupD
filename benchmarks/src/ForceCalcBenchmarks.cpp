/**
 * @file ForceCalcBenchmarks.cpp
 * @brief Benchmarks comparing different force calculation optimization strategies.
 *
 * Compares unoptimized O(n²), optimized (Newton's 3rd law), and alternative optimized
 * force calculation methods on SimpleContainer with LennardJonesForce.
 *
 * Filter: --benchmark_filter=Force/Calculation/
 */
#include <benchmark/benchmark.h>

#include "../code/simulationimpl/ForceCalcMethods.h"
#include "BenchmarkingUtils.h"
#include "particles/container/SimpleContainer.h"
#include "physics/LennardJonesForce.h"

namespace mol_sim {
/**
 * @brief Benchmarks the unoptimized force calculation (O(n²) full iteration).
 * Performs one calculateF() call on a container of 100-10,000 particles.
 */
template <ForceConcept forceType>
void bmForceCalcUnoptimized(benchmark::State& state) {
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
 * @brief Benchmarks the optimized force calculation (Newton's 3rd law).
 * Performs one calculateF() call on a container of 100-10,000 particles.
 */
template <ForceConcept forceType>
void bmForceCalcOptimized(benchmark::State& state) {
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
 * @brief Benchmarks the alternative optimized force calculation.
 * Performs one calculateF() call on a container of 100-10,000 particles.
 */
template <ForceConcept forceType>
void bmForceCalcOptimizedAlt(benchmark::State& state) {
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

BENCHMARK(bmForceCalcUnoptimized<LennardJonesForce>)
    ->Name("Force/Calculation/Unoptimized")
    ->RangeMultiplier(10)
    ->Range(100, 10000)
    ->Complexity()
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond)
    ->DisplayAggregatesOnly(true);
BENCHMARK(bmForceCalcOptimized<LennardJonesForce>)
    ->Name("Force/Calculation/Optimized")
    ->RangeMultiplier(10)
    ->Range(100, 10000)
    ->Complexity()
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond)
    ->DisplayAggregatesOnly(true);
BENCHMARK(bmForceCalcOptimizedAlt<LennardJonesForce>)
    ->Name("Force/Calculation/OptimizedAlt")
    ->RangeMultiplier(10)
    ->Range(100, 10000)
    ->Complexity()
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond)
    ->DisplayAggregatesOnly(true);
}  // namespace mol_sim
