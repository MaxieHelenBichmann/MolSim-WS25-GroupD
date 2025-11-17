#include <benchmark/benchmark.h>

#include <memory>

#include "../code/ForceCalcMethods.h"
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
BENCHMARK(bmForceOptimizedAlt<LennardJonesForce>)
    ->RangeMultiplier(10)
    ->Range(100, 1000)
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond);

/**
 * @brief Benchmarks the Simulation's calculateF method
 * Performs one calculateF() call on a container of 100-1000 particles using the Simulation class
 *
 */
void bmSimulationCalculateF(benchmark::State& state) {
    SimpleContainer particles;
    size_t n = state.range(0);
    particles.reserve(n);

    for (size_t i = 0; i < n; i++) {
        particles.addParticle(randomParticle());
    }

    SettingsParam settings(0.1, 0, 1, 5.0, 1.0);
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    Simulation<SimpleContainer> simulation(particles, std::move(force_source), settings, std::move(writer));

    for ([[maybe_unused]] auto _ : state) {
        simulation.calculateF();
        benchmark::DoNotOptimize(particles);
    }
}
BENCHMARK(bmSimulationCalculateF)
    ->RangeMultiplier(10)
    ->Range(100, 1000)
    ->Repetitions(10)
    ->Unit(benchmark::kMicrosecond);
}  // namespace mol_sim
