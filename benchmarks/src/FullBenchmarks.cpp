#include <benchmark/benchmark.h>

#include <random>

#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/ForceSource.h"
#include "utils/Simulation.h"
namespace mol_sim {
R3 randomR3() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> dist{};

    return {dist(gen), dist(gen), dist(gen)};
}

Particle randomParticle() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> dist{};

    return {randomR3(), randomR3(), dist(gen)};
}

void setUpLogging() { spdlog::set_level(spdlog::level::warn); }

static void BM_SimulationBig(benchmark::State& state) {
    setUpLogging();
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    size_t n = state.range(0);

    CuboidGenerator generator({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {std::move(n), std::move(n), std::move(n)}, 1.0, 1.0,
                              0.5);
    Simulation<SimpleContainer> simulation(part_container, GRAVITATIONAL, 0.1, 0, 1000);
    for (auto i : state) {
        generator.generateParticles(particles);
        simulation.run();
    }
    benchmark::DoNotOptimize(simulation);
    benchmark::DoNotOptimize(generator);
}
// Register the function as a benchmark
// BENCHMARK(BM_SimulationBig)->Range(4, 4 << 5)->Unit(benchmark::kMillisecond);

static void BM_SimulationGiven(benchmark::State& state) {
    setUpLogging();
    SimpleContainer part_container;
    ContainerRef particles(part_container);
    CuboidGenerator generator1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {40U, 8U, 1U}, 1.0, 1.0, 0.1);
    CuboidGenerator generator2({15.0, 15.0, 0.0}, {0.0, -10.0, 0.0}, {8U, 8U, 1U}, 1.0, 1.0, 0.1);
    Simulation<SimpleContainer> simulation(part_container, GRAVITATIONAL, 0.014, 0, 1000);
    for (auto i : state) {
        generator1.generateParticles(particles);
        generator2.generateParticles(particles);
        simulation.run();
    }
    benchmark::DoNotOptimize(simulation);
    benchmark::DoNotOptimize(generator1);
    benchmark::DoNotOptimize(generator2);
}
BENCHMARK(BM_SimulationGiven)->Unit(benchmark::kMillisecond);
;
BENCHMARK_MAIN();  //(NOLINT)
}  // namespace mol_sim
