/**
 * @file ContainerBenchmarks.cpp
 * @brief Benchmarks comparing concept-based (templated) vs abstract (virtual) container implementations.
 *
 * These benchmarks validate that compile-time polymorphism (C++ concepts) provides
 * better performance than runtime polymorphism (virtual functions) for particle containers.
 *
 * Filter: --benchmark_filter=Container/
 */
#include <benchmark/benchmark.h>

#include <memory>

#include "../code/containerimpl/AbstractContainer.h"
#include "../code/containerimpl/ContainerImpl.h"
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
void bmContainerConcept(benchmark::State& state) {
    containerType particles;
    size_t n = state.range(0);
    double res = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            particles.addParticle({static_cast<double>(i), 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 0., 0.);
        }
        benchmark::ClobberMemory();
        for (auto& p : particles) {
            res += p.getX()[0];
        }
        benchmark::DoNotOptimize(particles);
    }
    benchmark::DoNotOptimize(res);
};

/**
 * @brief Benchmarks the ContainerRef wrapper.
 * Runs the same benchmark as for the templated Container but wrapping the templated Container in a ContainerRef.
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
        benchmark::ClobberMemory();
        for (auto& p : particles) {
            res += p.getX()[0];
        }
        benchmark::DoNotOptimize(particles);
    }
    benchmark::DoNotOptimize(res);
};

/**
 * @brief Benchmarks a Container implemented using an abstract class (virtual functions).
 * Adds 8-8192 particles to the container, then iterates over all of them.
 * Reruns this Benchmark 10 times
 */
void bmContainerAbstract(benchmark::State& state) {
    std::unique_ptr<AbstractContainer> particles = std::make_unique<ContainerImpl>();
    size_t n = state.range(0);
    double res = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            particles->addParticle({static_cast<double>(i), 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 5., 1.);
        }
        benchmark::ClobberMemory();
        for (auto& p : *particles) {
            res += p.getX()[0];
        }
        benchmark::DoNotOptimize(particles);
    }
    benchmark::DoNotOptimize(res);
};

BENCHMARK(bmContainerConcept<SimpleContainer>)
    ->Name("Container/Concept/Small")
    ->Range(8 << 0, 8 << 6)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmContainerConcept<SimpleContainer>)
    ->Name("Container/Concept/Large")
    ->Range(16 << 6, 16 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(bmContainerRef<SimpleContainer>)
    ->Name("Container/Ref/Small")
    ->Range(8 << 0, 8 << 6)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmContainerRef<SimpleContainer>)
    ->Name("Container/Ref/Large")
    ->Range(16 << 6, 16 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(bmContainerAbstract)
    ->Name("Container/Abstract/Small")
    ->Range(8 << 0, 8 << 6)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmContainerAbstract)
    ->Name("Container/Abstract/Large")
    ->Range(16 << 6, 16 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

}  // namespace mol_sim
