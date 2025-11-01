#include <benchmark/benchmark.h>

#include "particles/Particle.h"

static void BM_ParticleCreation(benchmark::State& state) {
    for (auto _ : state) {
    }
}
BENCHMARK(BM_ParticleCreation);

BENCHMARK_MAIN();
