#include <benchmark/benchmark.h>

#include <memory>

#include "particles/Particle.h"
#include "physics/ForceSource.h"
#include "physics/GravitationalForce.h"

namespace mol_sim {
/**
 * @brief Benchmarks a templated Container of Type ContainerType.
 * Adds 8-8192 particles to the container, then iterates over all off them.
 * Reruns this Benchmark 10 Times
 * @tparam containerType Type of container to be benchmarked
 */

void bmAbstractForce(benchmark::State& state) {
    std::unique_ptr<ForceSource> force = std::make_unique<GravitationalForce>();
    Particle p1 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};
    Particle p2 = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            [[maybe_unused]] auto x = force->calculateForce(p1, p2);
        }
    }
};

BENCHMARK(bmAbstractForce)
    ->RangeMultiplier(10)
    ->Range(1000, 1000000)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

}  // namespace mol_sim
