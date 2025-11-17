#include <benchmark/benchmark.h>

#include <memory>

#include "../code/AbstractForce.h"
#include "../code/ForceConcept.h"
#include "../code/GravitationalAbstract.h"
#include "../code/GravitationalConcept.h"
#include "particles/Particle.h"
#include "physics/GravitationalForce.h"

namespace mol_sim {
/**
 * @brief Benchmarks a force Source implemented using an abstract class.
 * Runs the gravitational force calculation 1.000-1.000.000 times
 * Reruns this Benchmark 10 Times
 */

void bmAbstractForce(benchmark::State& state) {
    std::unique_ptr<ForceAbstract> force = std::make_unique<GravitationalAbstract>();
    Particle p1 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    Particle p2 = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            [[maybe_unused]] auto x = force->applyForce(p1, p2);
        }
    }
};

BENCHMARK(bmAbstractForce)
    ->RangeMultiplier(10)
    ->Range(1000, 1000000)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
/**
 * @brief Benchmarks a force Source implemented using a concept.
 * Runs the gravitational force calculation 1.000-1.000.000 times
 * Reruns this Benchmark 10 Times
 * @tparam forceType Type of force to benchmark
 */
template <ForceConcept forceType>
void bmConceptForce(benchmark::State& state) {
    forceType force;
    Particle p1 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    Particle p2 = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 5., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            [[maybe_unused]] auto x = force.applyForce(p1, p2);
        }
    }
};

BENCHMARK(bmConceptForce<GravitationalConcept>)
    ->RangeMultiplier(10)
    ->Range(1000, 1000000)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

}  // namespace mol_sim
