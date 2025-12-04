/**
 * @file CellBenchmark.cpp
 * @brief Benchmarks comparing different cell data structure implementations.
 *
 * Compares Vector, Set, and UnorderedSet implementations for cell operations:
 * add, iterate, update, remove, and mixed operations. Used to determine the
 * optimal data structure for LinkedCellContainer's internal cells.
 *
 * Filter: --benchmark_filter=Cell/
 */
#include <benchmark/benchmark.h>

#include <array>

#include "../code/cellimpl/CellSet.h"
#include "../code/cellimpl/CellUnorderedSet.h"
#include "../code/cellimpl/CellVector.h"
namespace mol_sim {

/**
 * @brief Benchmarks particle addition using std::vector-backed cell.
 * Adds N particles to a cell (N = 8 to 8192).
 */
void bmCellAddVector(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle addition using std::set-backed cell.
 * Adds N particles to a cell (N = 8 to 8192).
 */
void bmCellAddSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle addition using std::unordered_set-backed cell.
 * Adds N particles to a cell (N = 8 to 8192).
 */
void bmCellAddUnorderedSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks iteration over std::vector-backed cell.
 * Iterates over N particles in a cell (N = 8 to 8192).
 */
void bmCellIteratorVector(benchmark::State& state) {
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

/**
 * @brief Benchmarks iteration over std::set-backed cell.
 * Iterates over N particles in a cell (N = 8 to 8192).
 */
void bmCellIteratorSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks iteration over std::unordered_set-backed cell.
 * Iterates over N particles in a cell (N = 8 to 8192).
 */
void bmCellIteratorUnorderedSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle index updates in std::vector-backed cell.
 * Updates N particle indices (N = 8 to 8192).
 */
void bmCellUpdateVector(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle index updates in std::set-backed cell.
 * Updates N particle indices (N = 8 to 8192).
 */
void bmCellUpdateSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle index updates in std::unordered_set-backed cell.
 * Updates N particle indices (N = 8 to 8192).
 */
void bmCellUpdateUnorderedSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle removal from std::vector-backed cell.
 * Removes N particles from a cell (N = 8 to 8192).
 */
void bmCellRemoveVector(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle removal from std::set-backed cell.
 * Removes N particles from a cell (N = 8 to 8192).
 */
void bmCellRemoveSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks particle removal from std::unordered_set-backed cell.
 * Removes N particles from a cell (N = 8 to 8192).
 */
void bmCellRemoveUnorderedSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks mixed operations on std::vector-backed cell.
 * Performs add, iterate, and remove operations (N = 8 to 8192).
 */
void bmCellMixedVector(benchmark::State& state) {
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

/**
 * @brief Benchmarks mixed operations on std::set-backed cell.
 * Performs add, iterate, and remove operations (N = 8 to 8192).
 */
void bmCellMixedSet(benchmark::State& state) {
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

/**
 * @brief Benchmarks mixed operations on std::unordered_set-backed cell.
 * Performs add, iterate, and remove operations (N = 8 to 8192).
 */
void bmCellMixedUnorderedSet(benchmark::State& state) {
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

BENCHMARK(bmCellAddVector)
    ->Name("Cell/Add/Vector")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellAddSet)
    ->Name("Cell/Add/Set")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellAddUnorderedSet)
    ->Name("Cell/Add/UnorderedSet")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmCellIteratorVector)
    ->Name("Cell/Iterator/Vector")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellIteratorSet)
    ->Name("Cell/Iterator/Set")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellIteratorUnorderedSet)
    ->Name("Cell/Iterator/UnorderedSet")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmCellUpdateVector)
    ->Name("Cell/Update/Vector")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellUpdateSet)
    ->Name("Cell/Update/Set")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellUpdateUnorderedSet)
    ->Name("Cell/Update/UnorderedSet")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmCellRemoveVector)
    ->Name("Cell/Remove/Vector")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellRemoveSet)
    ->Name("Cell/Remove/Set")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
BENCHMARK(bmCellRemoveUnorderedSet)
    ->Name("Cell/Remove/UnorderedSet")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(bmCellMixedVector)
    ->Name("Cell/Mixed/Vector")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmCellMixedSet)
    ->Name("Cell/Mixed/Set")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(bmCellMixedUnorderedSet)
    ->Name("Cell/Mixed/UnorderedSet")
    ->RangeMultiplier(2)
    ->Range(8, 8 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);

}  // namespace mol_sim
