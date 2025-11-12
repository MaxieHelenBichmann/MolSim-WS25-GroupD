#include <benchmark/benchmark.h>

#include "../code/ForceCalcMethods.h"
#include "BenchmarkingUtils.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "physics/LennardJonesForce.h"

namespace mol_sim {
/**
 * @brief Benchmarks the unoptimized Force Calculation
 * Performs one calulateF() call on a container of 100-1000 particles
 *
 */
template <ForceSource forceType>
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
BENCHMARK(bmForceUnOptimized<LennardJonesForce>)
    ->RangeMultiplier(10)
    ->Range(100, 1000)
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond);
/**
 * @brief Benchmarks the optimized Force Calculation
 * Performs one calulateF() call on a container of 100-1000 particles
 *
 */
template <ForceSource forceType>
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
BENCHMARK(bmForceOptimized<LennardJonesForce>)
    ->RangeMultiplier(10)
    ->Range(100, 1000)
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond);

/**
 * @brief Benchmarks the alternative optimized Force Calculation
 * Performs one calulateF() call on a container of 100-1000 particles
 *
 */
template <ForceSource forceType>
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
BENCHMARK(bmForceOptimizedAlt<LennardJonesForce>)
    ->RangeMultiplier(10)
    ->Range(100, 1000)
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond);
}  // namespace mol_sim
