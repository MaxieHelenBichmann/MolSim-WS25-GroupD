#include <benchmark/benchmark.h>

#include <memory>

#include "../code/ForceCalcMethods.h"
#include "BenchmarkingUtils.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "physics/GravitationalForce.h"

namespace mol_sim {
/**
 * @brief Benchmarks the unoptimized Force Calculation
 *
 *
 */
void bmForceUnOptimized(benchmark::State& state) {
    SimpleContainer particles;
    size_t n = state.range(0);
    particles.reserve(n);
    std::unique_ptr<ForceSource> force_source = std::make_unique<GravitationalForce>();
    for (size_t i = 0; i < n; i++) {
        particles.addParticle(randomParticle());
    }
    for ([[maybe_unused]] auto _ : state) {
        calculateFUnoptimized<SimpleContainer>(particles, force_source);
        benchmark::DoNotOptimize(particles);
    }
}
BENCHMARK(bmForceUnOptimized)->RangeMultiplier(10)->Range(100, 1000)->Repetitions(10)->Unit(benchmark::kMicrosecond);
/**
 * @brief Benchmarks the optimized Force Calculation
 *
 *
 */
void bmForceOptimized(benchmark::State& state) {
    SimpleContainer particles;
    size_t n = state.range(0);
    particles.reserve(n);
    std::unique_ptr<ForceSource> force_source = std::make_unique<GravitationalForce>();
    for (size_t i = 0; i < n; i++) {
        particles.addParticle(randomParticle());
    }

    for ([[maybe_unused]] auto _ : state) {
        calculateF<SimpleContainer>(particles, force_source);
        benchmark::DoNotOptimize(particles);
    }
}
BENCHMARK(bmForceOptimized)->RangeMultiplier(10)->Range(100, 1000)->Repetitions(10)->Unit(benchmark::kMicrosecond);
}  // namespace mol_sim
