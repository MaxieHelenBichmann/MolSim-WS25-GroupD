#include <benchmark/benchmark.h>

#include <memory>

#include "../code/AbstractContainer.h"
#include "../code/ContainerImpl.h"
#include "BenchmarkingUtils.h"
#include "particles/ParticleContainer.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"

namespace mol_sim {
/**
 * @brief Benchmarks a templated Container of Type ContainerType.
 * Adds 8-8192 particles to the container, then iterates over all of them.
 * Reruns this Benchmark 10 times
 * @tparam containerType Type of container to be benchmarked
 */
template <ParticleContainer containerType>
void bmTemplatedContainer(benchmark::State& state) {
    containerType particles;
    size_t n = state.range(0);
    double res = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            particles.addParticle({static_cast<double>(i), 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 0., 0.);
        }
        for (auto& p : particles) {
            res += p.getX()[0];
        }
    }
    benchmark::DoNotOptimize(res);
};
BENCHMARK(bmTemplatedContainer<SimpleContainer>)
    ->Range(8 << 0, 8 << 6)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmTemplatedContainer<SimpleContainer>)
    ->Range(16 << 6, 16 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
;
/**
 * @brief Benchmarks the ContainerRef.
 * By running the same Benchmark as for the templated Container but wrapping the templated Container in a ContainerRef.
 * @tparam containerType Type of container to be wrapped
 */
template <ParticleContainer containerType>
void bmContainerRef(benchmark::State& state) {
    containerType part_container;
    ContainerRef particles(part_container);
    size_t n = state.range(0);
    double res = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            particles.addParticle({static_cast<double>(i), 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 0., 0.);
        }
        for (auto& p : particles) {
            res += p.getX()[0];
        }
    }
    benchmark::DoNotOptimize(res);
};

BENCHMARK(bmContainerRef<SimpleContainer>)
    ->Range(8 << 0, 8 << 6)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmContainerRef<SimpleContainer>)
    ->Range(16 << 6, 16 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

void bmAbstractContainer(benchmark::State& state) {
    std::unique_ptr<AbstractContainer> particles = std::make_unique<ContainerImpl>();
    size_t n = state.range(0);
    double res = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            particles->addParticle({static_cast<double>(i), 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 5., 1.);
        }
        for (auto& p : *particles) {
            res += p.getX()[0];
        }
    }
    benchmark::DoNotOptimize(res);
};
BENCHMARK(bmAbstractContainer)
    ->Range(8 << 0, 8 << 6)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmAbstractContainer)
    ->Range(16 << 6, 16 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

}  // namespace mol_sim
