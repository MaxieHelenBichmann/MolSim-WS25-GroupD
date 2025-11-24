#include <benchmark/benchmark.h>

#include <array>

#include "../code/CellSet.h"
#include "../code/CellUnorderedSet.h"
#include "../code/CellVector.h"
namespace mol_sim {

// Add particle Benchmark

void bmAddVector(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    CellVector cell(bounds);
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
    }
}
void bmAddSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    CellSet cell(bounds);
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
    }
}
void bmAddUnOrdSet(benchmark::State& state) {
    std::array<double, 6> bounds = {0., 1., 0., 1., 0., 1.};
    CellUnorderedSet cell(bounds);
    size_t n = state.range(0);
    for ([[maybe_unused]] auto _ : state) {
        for (size_t i = 0; i < n; i++) {
            cell.addParticle(i);
        }
    }
}

// Iterator Benchmark

// Update Benchmark

// Combined
}  // namespace mol_sim
