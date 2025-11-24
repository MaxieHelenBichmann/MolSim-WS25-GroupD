#include <benchmark/benchmark.h>

#include <array>

#include "../code/cellimpl/CellSet.h"
#include "../code/cellimpl/CellUnorderedSet.h"
#include "../code/cellimpl/CellVector.h"
namespace mol_sim {

// Add particle Benchmark

void bmAddVector(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        CellVector cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        benchmark::DoNotOptimize(cell);
    }
}
void bmAddSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        CellSet cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        benchmark::DoNotOptimize(cell);
    }
}
void bmAddUnOrdSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        CellUnorderedSet cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        benchmark::DoNotOptimize(cell);
    }
}

// Iterator Benchmark

void bmIteratorVector(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    CellVector cell(bounds);
    size_t n = state.range(0);

    for (size_t i = 0; i < n; i++) {
        cell.addParticle(i);
    }
    size_t res = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (auto p : cell.particles()) {
            res += p;
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(res);
}

void bmIteratorSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    CellSet cell(bounds);
    size_t n = state.range(0);

    for (size_t i = 0; i < n; i++) {
        cell.addParticle(i);
    }
    size_t res = 0;
    for ([[maybe_unused]] auto _ : state) {
        for (auto p : cell.particles()) {
            res += p;
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(res);
}

void bmIteratorUnOrdSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    CellUnorderedSet cell(bounds);
    size_t n = state.range(0);

    for (size_t i = 0; i < n; i++) {
        cell.addParticle(i);
    }
    size_t res = 0;

    for ([[maybe_unused]] auto _ : state) {
        auto it = cell.stableIteratorBegin();
        auto end = cell.stableIteratorEnd();
        for (; it != end; ++it) {
            res += *it;
        }
        benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(res);
}

// Update Benchmark

void bmUpdateVector(benchmark::State& state) {
    size_t n = state.range(0);
    std::array<double, 6> bounds = {-static_cast<double>(n), static_cast<double>(n),  -static_cast<double>(n),
                                    static_cast<double>(n),  -static_cast<double>(n), static_cast<double>(n)};
    CellVector cell(bounds);
    for (size_t i = 0; i < n; i++) {
        cell.addParticle(i);
    }

    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            cell.updateParticleIndex(i, -i);
        }
        benchmark::DoNotOptimize(cell);
    }
}

void bmUpdateSet(benchmark::State& state) {
    size_t n = state.range(0);
    std::array<double, 6> bounds = {-static_cast<double>(n), static_cast<double>(n),  -static_cast<double>(n),
                                    static_cast<double>(n),  -static_cast<double>(n), static_cast<double>(n)};
    CellSet cell(bounds);
    for (size_t i = 0; i < n; i++) {
        cell.addParticle(i);
    }

    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            cell.updateParticleIndex(i, -i);
        }
        benchmark::DoNotOptimize(cell);
    }
}

void bmUpdateUnOrdSet(benchmark::State& state) {
    size_t n = state.range(0);
    std::array<double, 6> bounds = {-static_cast<double>(n), static_cast<double>(n),  -static_cast<double>(n),
                                    static_cast<double>(n),  -static_cast<double>(n), static_cast<double>(n)};
    CellUnorderedSet cell(bounds);
    for (size_t i = 0; i < n; i++) {
        cell.addParticle(i);
    }
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            cell.updateParticleIndex(i, -i);
        }
        benchmark::DoNotOptimize(cell);
    }
}

// Remove Benchmark

void bmRemoveVector(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        CellVector cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        state.ResumeTiming();
        for (size_t i = 0; i < n; i++) {
            cell.removeParticle(i);
        }
    }
}

void bmRemoveSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        CellSet cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        state.ResumeTiming();
        for (size_t i = 0; i < n; i++) {
            cell.removeParticle(i);
        }
    }
}

void bmRemoveUnOrdSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        CellUnorderedSet cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        state.ResumeTiming();
        for (size_t i = 0; i < n; i++) {
            cell.removeParticle(i);
        }
    }
}

// Mixed Operations Benchmark

void bmMixedVector(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);

    for ([[maybe_unused]] auto _ : state) {
        CellVector cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        size_t sum = 0;
        for (auto p : cell.particles()) {
            sum += p;
        }
        for (size_t i = 0; i < n / 2; i++) {
            cell.removeParticle(i);
        }
        benchmark::DoNotOptimize(cell);
        benchmark::DoNotOptimize(sum);
    }
}

void bmMixedSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);

    for ([[maybe_unused]] auto _ : state) {
        CellSet cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        size_t sum = 0;
        for (auto p : cell.particles()) {
            sum += p;
        }
        for (size_t i = 0; i < n / 2; i++) {
            cell.removeParticle(i);
        }
        benchmark::DoNotOptimize(cell);
        benchmark::DoNotOptimize(sum);
    }
}

void bmMixedUnOrdSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    size_t n = state.range(0);

    for ([[maybe_unused]] auto _ : state) {
        CellUnorderedSet cell(bounds);
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
        size_t sum = 0;
        auto it = cell.stableIteratorBegin();
        auto end = cell.stableIteratorEnd();
        for (; it != end; ++it) {
            sum += *it;
        }
        for (size_t i = 0; i < n / 2; i++) {
            cell.removeParticle(i);
        }
        benchmark::DoNotOptimize(cell);
        benchmark::DoNotOptimize(sum);
    }
}

// Benchmark Registrations

BENCHMARK(bmAddVector)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmAddSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmAddUnOrdSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmIteratorVector)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmIteratorSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmIteratorUnOrdSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmUpdateVector)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmUpdateSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmUpdateUnOrdSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmRemoveVector)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmRemoveSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmRemoveUnOrdSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmMixedVector)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmMixedSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmMixedUnOrdSet)
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);

}  // namespace mol_sim
