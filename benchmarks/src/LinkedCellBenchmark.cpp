#include <benchmark/benchmark.h>

#include <cstdint>

#include "../code/linkedcellimpl/LinkedCellContainerDirect.h"
#include "../code/linkedcellimpl/LinkedCellContainerExplicit.h"
#include "../src/BenchmarkingUtils.h"

namespace mol_sim {

// Add particle Benchmark

void bmAddDirect(benchmark::State& state) {
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

void bmAddExplicit(benchmark::State& state) {
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

// Iterator Benchmark (all particles)

void bmIteratorDirect(benchmark::State& state) {
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

void bmIteratorExplicit(benchmark::State& state) {
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

// Proximity Iterator Benchmark

void bmProximityIteratorDirect(benchmark::State& state) {
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

    R3 center = {50.0, 50.0, 50.0};
    double radius = 10.0;
    size_t count = 0;

    for ([[maybe_unused]] auto _ : state) {
        auto it = container.proximityBegin(center, radius);
        auto end = container.proximityEnd(center, radius);
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

void bmProximityIteratorExplicit(benchmark::State& state) {
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

    R3 center = {50.0, 50.0, 50.0};
    double radius = 10.0;
    size_t count = 0;

    for ([[maybe_unused]] auto _ : state) {
        auto it = container.proximityBegin(center, radius);
        auto end = container.proximityEnd(center, radius);
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

// Update Position Benchmark (with cell switching)

void bmUpdatePositionDirect(benchmark::State& state) {
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
                container.updateParticlePosition(*it, new_pos);
            }
        }
        benchmark::DoNotOptimize(container);
    }
}

void bmUpdatePositionExplicit(benchmark::State& state) {
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

// Clear Benchmark

void bmClearDirect(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        LinkedCellContainerDirect container(domain_size, cutoff);
        for (size_t i = 0; i < n; i++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }
        state.ResumeTiming();
        container.clear();
        benchmark::DoNotOptimize(container);
    }
}

void bmClearExplicit(benchmark::State& state) {
    size_t n = state.range(0);
    R3 domain_size = {100.0, 100.0, 100.0};
    double cutoff = 3.0;
    std::array<double, 6> bounds = {0., domain_size[0], 0., domain_size[1], 0., domain_size[2]};

    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        LinkedCellContainerExplicit container(domain_size, cutoff);
        for (size_t i = 0; i < n; i++) {
            R3 pos = randomR3(bounds);
            R3 vel = randomR3();
            container.addParticle(pos, vel, 1.0, 1.0, 1.0);
        }
        state.ResumeTiming();
        container.clear();
        benchmark::DoNotOptimize(container);
    }
}

// Mixed Operations Benchmark

void bmMixedDirect(benchmark::State& state) {
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

        size_t count = 0;
        for (auto& p : container) {
            count += p.getType();
        }

        R3 center = {50.0, 50.0, 50.0};
        double radius = 10.0;
        auto it = container.proximityBegin(center, radius);
        auto end = container.proximityEnd(center, radius);
        for (; it != end; ++it) {
            count += it->getType();
        }

        benchmark::DoNotOptimize(container);
        benchmark::DoNotOptimize(count);
    }
}

void bmMixedExplicit(benchmark::State& state) {
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

        size_t count = 0;
        for (auto& p : container) {
            count += p.getType();
        }

        R3 center = {50.0, 50.0, 50.0};
        double radius = 10.0;
        auto it = container.proximityBegin(center, radius);
        auto end = container.proximityEnd(center, radius);
        for (; it != end; ++it) {
            count += it->getType();
        }

        benchmark::DoNotOptimize(container);
        benchmark::DoNotOptimize(count);
    }
}

// Boundary Iterator Benchmark

void bmBoundaryIteratorDirect(benchmark::State& state) {
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

void bmBoundaryIteratorExplicit(benchmark::State& state) {
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

// Halo Iterator Benchmark

void bmHaloIteratorDirect(benchmark::State& state) {
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
        auto it = container.haloBegin();
        auto end = container.haloEnd();
        for (; it != end; ++it) {
            count += it->getType();
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(count);
}

void bmHaloIteratorExplicit(benchmark::State& state) {
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

BENCHMARK(bmAddDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();
BENCHMARK(bmAddExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond)
    ->Complexity();

BENCHMARK(bmIteratorDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();
BENCHMARK(bmIteratorExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();

BENCHMARK(bmProximityIteratorDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();
BENCHMARK(bmProximityIteratorExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();

BENCHMARK(bmUpdatePositionDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(bmUpdatePositionExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(bmClearDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();
BENCHMARK(bmClearExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond)
    ->Complexity();

BENCHMARK(bmBoundaryIteratorDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(bmBoundaryIteratorExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(bmHaloIteratorDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(bmHaloIteratorExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(bmMixedDirect)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(bmMixedExplicit)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMillisecond);

}  // namespace mol_sim
