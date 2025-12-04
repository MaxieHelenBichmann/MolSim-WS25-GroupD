/**
 * @file ForceBenchmark.cpp
 * @brief Benchmarks comparing concept-based vs abstract (virtual) force implementations.
 *
 * These benchmarks validate that compile-time polymorphism (C++ concepts) provides
 * better performance than runtime polymorphism (virtual functions) for force calculations.
 *
 * Filter: --benchmark_filter=Force/Polymorphism/
 */
#include <benchmark/benchmark.h>

#include <memory>

#include "../code/forceimpl/AbstractForce.h"
#include "../code/forceimpl/ForceConcept.h"
#include "../code/forceimpl/GravitationalAbstract.h"
#include "../code/forceimpl/GravitationalConcept.h"
#include "particles/Particle.h"

namespace mol_sim {
/**
 * @brief Benchmarks a force source implemented using an abstract class (virtual functions).
 * Runs the gravitational force calculation 1,000-1,000,000 times.
 * Reruns this Benchmark 10 times.
 */
void bmForceAbstract(benchmark::State& state) {
    std::unique_ptr<ForceAbstract> force = std::make_unique<GravitationalAbstract>();
    Particle p1 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    Particle p2 = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    benchmark::DoNotOptimize(p1);
    benchmark::DoNotOptimize(p2);
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            auto x = force->applyForce(p1, p2);
            benchmark::DoNotOptimize(x);
        }
    }
};

/**
 * @brief Benchmarks a force source implemented using a concept (compile-time polymorphism).
 * Runs the gravitational force calculation 1,000-1,000,000 times.
 * Reruns this Benchmark 10 times.
 * @tparam forceType Type of force to benchmark
 */
template <ForceConcept forceType>
void bmForceConcept(benchmark::State& state) {
    forceType force;
    Particle p1 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    Particle p2 = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    benchmark::DoNotOptimize(p1);
    benchmark::DoNotOptimize(p2);
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            auto x = force.applyForce(p1, p2);
            benchmark::DoNotOptimize(x);
        }
    }
};

BENCHMARK(bmForceAbstract)
    ->Name("Force/Polymorphism/Abstract")
    ->RangeMultiplier(10)
    ->Range(1000, 1000000)
    ->Complexity()
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(bmForceConcept<GravitationalConcept>)
    ->Name("Force/Polymorphism/Concept")
    ->RangeMultiplier(10)
    ->Range(1000, 1000000)
    ->Complexity()
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
}  // namespace mol_sim
