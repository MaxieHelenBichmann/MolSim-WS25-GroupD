#include <benchmark/benchmark.h>

namespace mol_sim {
/**
 * @brief Benchmarks a templated Container of Type ContainerType.
 * Adds 8-8192 particles to the container, then iterates over all off them.
 * Reruns this Benchmark 10 Times
 * @tparam containerType Type of container to be benchmarked
 */

void BM_AbstractForce(benchmark::State& state) {

};
/*
BENCHMARK(BM_AbstractForce)
    ->Range(8 << 0, 8 << 6)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_AbstractForce)
    ->Range(16 << 6, 16 << 10)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true)
    ->Unit(benchmark::kMicrosecond);
*/

}  // namespace mol_sim
