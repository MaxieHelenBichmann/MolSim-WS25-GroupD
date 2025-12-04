/**
 * @file LinkedCellBenchmark.cpp
 * @brief Benchmarks comparing LinkedCellContainer implementation strategies.
 *
 * Compares Direct (iterator-based cell management) vs Explicit (index-based cell management)
 * implementations across various operations: add, iterate, proximity, update, boundary, halo.
 *
 * Filter: --benchmark_filter=LinkedCell/
 */
#include <benchmark/benchmark.h>

#include "../code/linkedcellimpl/LinkedCellContainerDirect.h"
#include "../code/linkedcellimpl/LinkedCellContainerExplicit.h"
#include "../src/BenchmarkingUtils.h"

namespace mol_sim {

/**
 * @brief Benchmarks particle addition to LinkedCellContainerDirect.
 * Adds N random particles (N = 128 to 8192).
 */
void bmLinkedCellAddDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    for ([[maybe_unused]] auto _ : state) {
        LinkedCellContainerDirect container(domain_size, cutoff);
        for (size_t i = 0; i < n; i++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }
        benchmark::DoNotOptimize(container);
    }
}

/**
 * @brief Benchmarks particle addition to LinkedCellContainerExplicit.
 * Adds N random particles (N = 128 to 8192).
 */
void bmLinkedCellAddExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    for ([[maybe_unused]] auto _ : state) {
        LinkedCellContainerExplicit container(domain_size, cutoff);
        for (size_t i = 0; i < n; i++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }
        benchmark::DoNotOptimize(container);
    }
}

/**
 * @brief Benchmarks full iteration over LinkedCellContainerDirect.
 * Iterates over N particles (N = 128 to 8192).
 */
void bmLinkedCellIteratorDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerDirect container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    size_t count = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (auto& p : container) {
            count += p.getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

/**
 * @brief Benchmarks full iteration over LinkedCellContainerExplicit.
 * Iterates over N particles (N = 128 to 8192).
 */
void bmLinkedCellIteratorExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerExplicit container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    size_t count = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (auto& p : container) {
            count += p.getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

/**
 * @brief Benchmarks proximity iteration on LinkedCellContainerDirect.
 * Finds neighbors within cutoff radius from center position (N = 128 to 8192).
 */
void bmLinkedCellProximityDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 10.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerDirect container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    R3 center = {50.0, 50.0, 50.0};
    size_t count = 0;

    for ([[maybe_unused]] auto _ : state) {
        auto it = container.proximityBegin(center);
        auto end = container.proximityEnd(center);
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

/**
 * @brief Benchmarks proximity iteration on LinkedCellContainerExplicit.
 * Finds neighbors within cutoff radius from center position (N = 128 to 8192).
 */
void bmLinkedCellProximityExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 10.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerExplicit container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    R3 center = {50.0, 50.0, 50.0};
    size_t count = 0;

    for ([[maybe_unused]] auto _ : state) {
        auto it = container.proximityBegin(center);
        auto end = container.proximityEnd(center);
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

/**
 * @brief Benchmarks position updates on LinkedCellContainerDirect.
 * Updates particle positions causing cell migrations (N = 128 to 8192).
 */
void bmLinkedCellUpdateDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerDirect container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    for ([[maybe_unused]] auto _ : state) {
        for (auto it = container.begin(); it != container.end(); ++it) {
            R3 new_pos = it->getX() + R3{1.5, 1.5, 1.5};
            if (container.fitsContainer(new_pos)) {
                container.updateParticlePosition(it, new_pos);
            }
        }
        benchmark::DoNotOptimize(container);
    }
}

/**
 * @brief Benchmarks position updates on LinkedCellContainerExplicit.
 * Updates particle positions causing cell migrations (N = 128 to 8192).
 */
void bmLinkedCellUpdateExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerExplicit container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    for ([[maybe_unused]] auto _ : state) {
        for (auto it = container.begin(); it != container.end(); ++it) {
            R3 new_pos = it->getX() + R3{1.5, 1.5, 1.5};
            if (container.fitsContainer(new_pos)) {
                container.updateParticlePosition(it, new_pos);
            }
        }
        benchmark::DoNotOptimize(container);
    }
}

/**
 * @brief Benchmarks mixed operations on LinkedCellContainerDirect.
 * Performs add, iterate, and proximity operations (N = 128 to 8192).
 */
void bmLinkedCellMixedDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 10.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    for ([[maybe_unused]] auto _ : state) {
        LinkedCellContainerDirect container(domain_size, cutoff);
        for (size_t i = 0; i < n; i++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }

        size_t count = 0;
        for (auto& p : container) {
            count += p.getType();
        }

        R3 center = {50.0, 50.0, 50.0};
        auto it = container.proximityBegin(center);
        auto end = container.proximityEnd(center);
        for (; it != end; ++it) {
            count += it->getType();
        }

        benchmark::DoNotOptimize(container);
        benchmark::DoNotOptimize(count);
    }
}

/**
 * @brief Benchmarks mixed operations on LinkedCellContainerExplicit.
 * Performs add, iterate, and proximity operations (N = 128 to 8192).
 */
void bmLinkedCellMixedExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 10.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    for ([[maybe_unused]] auto _ : state) {
        LinkedCellContainerExplicit container(domain_size, cutoff);
        for (size_t i = 0; i < n; i++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }

        size_t count = 0;
        for (auto& p : container) {
            count += p.getType();
        }

        R3 center = {50.0, 50.0, 50.0};
        auto it = container.proximityBegin(center);
        auto end = container.proximityEnd(center);
        for (; it != end; ++it) {
            count += it->getType();
        }

        benchmark::DoNotOptimize(container);
        benchmark::DoNotOptimize(count);
    }
}

/**
 * @brief Benchmarks boundary iteration on LinkedCellContainerDirect.
 * Iterates over particles in boundary cells (N = 128 to 8192).
 */
void bmLinkedCellBoundaryDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerDirect container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    size_t count = 0;
    for ([[maybe_unused]] auto _ : state) {
        auto it = container.boundaryBegin();
        auto end = container.boundaryEnd();
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

/**
 * @brief Benchmarks boundary iteration on LinkedCellContainerExplicit.
 * Iterates over particles in boundary cells (N = 128 to 8192).
 */
void bmLinkedCellBoundaryExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerExplicit container(domain_size, cutoff);
    for (size_t i = 0; i < n; i++) {
        R3 pos = randomR3(bounds);
        R3 vel = randomR3();
        container.addParticle(pos, vel, 1.0, 1.0, 1.0);
    }

    size_t count = 0;
    for ([[maybe_unused]] auto _ : state) {
        auto it = container.boundaryBegin();
        auto end = container.boundaryEnd();
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

/**
 * @brief Benchmarks halo iteration on LinkedCellContainerDirect.
 * Iterates over particles in halo cells (N = 128 to 8192).
 */
void bmLinkedCellHaloDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};
    LinkedCellContainerDirect container(domain_size, cutoff);

    for (size_t i = 0; i < 6; i++) {
        bounds[i] += (2.0 * static_cast<double>(i % 2)) - 1.0;
        for (size_t j = 0; j < n / 6; j++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }
    }

    size_t count = 0;
    for ([[maybe_unused]] auto _ : state) {
        auto it = container.haloBegin();
        auto end = container.haloEnd();
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

/**
 * @brief Benchmarks halo iteration on LinkedCellContainerExplicit.
 * Iterates over particles in halo cells (N = 128 to 8192).
 */
void bmLinkedCellHaloExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    LinkedCellContainerExplicit container(domain_size, cutoff);
    for (size_t i = 0; i < 6; i++) {
        bounds[i] += (2.0 * static_cast<double>(i % 2)) - 1.0;
        for (size_t j = 0; j < n / 6; j++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }
    }

    size_t count = 0;
    for ([[maybe_unused]] auto _ : state) {
        auto it = container.haloBegin();
        auto end = container.haloEnd();
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

// Benchmark Registrations

BENCHMARK(bmLinkedCellAddDirect)
    ->Name("LinkedCell/Add/Direct")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();
BENCHMARK(bmLinkedCellAddExplicit)
    ->Name("LinkedCell/Add/Explicit")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();

BENCHMARK(bmLinkedCellIteratorDirect)
    ->Name("LinkedCell/Iterator/Direct")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();
BENCHMARK(bmLinkedCellIteratorExplicit)
    ->Name("LinkedCell/Iterator/Explicit")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();

// there seems to be be a problem with the direct proximity iterators immediately returning
BENCHMARK(bmLinkedCellProximityDirect)
    ->Name("LinkedCell/Proximity/Direct")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();
// there seems to be be a problem with the explicit proximity iterators causing segfaults
BENCHMARK(bmLinkedCellProximityExplicit)
    ->Name("LinkedCell/Proximity/Explicit")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();

BENCHMARK(bmLinkedCellUpdateDirect)
    ->Name("LinkedCell/Update/Direct")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(bmLinkedCellUpdateExplicit)
    ->Name("LinkedCell/Update/Explicit")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(bmLinkedCellBoundaryDirect)
    ->Name("LinkedCell/Boundary/Direct")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(bmLinkedCellBoundaryExplicit)
    ->Name("LinkedCell/Boundary/Explicit")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(bmLinkedCellHaloDirect)
    ->Name("LinkedCell/Halo/Direct")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(bmLinkedCellHaloExplicit)
    ->Name("LinkedCell/Halo/Explicit")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(bmLinkedCellMixedDirect)
    ->Name("LinkedCell/Mixed/Direct")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(bmLinkedCellMixedExplicit)
    ->Name("LinkedCell/Mixed/Explicit")
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);

}  // namespace mol_sim
