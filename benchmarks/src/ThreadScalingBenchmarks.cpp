/**
 * @file ThreadScalingBenchmarks.cpp
 * @brief Benchmarks for OpenMP thread scaling performance.
 *
 * Measures parallel efficiency across different thread counts for realistic workloads.
 * Tests strong scaling (fixed problem size) and weak scaling (problem size grows with threads).
 *
 * Usage:
 *   # Run all thread scaling benchmarks
 *   ./build/benchmarks/MolSimBench --benchmark_filter=ThreadScaling
 *
 *   # Run only strong scaling (fixed particle count)
 *   ./build/benchmarks/MolSimBench --benchmark_filter=ThreadScaling/Strong
 *
 *   # Run only weak scaling (particles scale with threads)
 *   ./build/benchmarks/MolSimBench --benchmark_filter=ThreadScaling/Weak
 *
 *   # Export results to JSON for analysis
 *   ./build/benchmarks/MolSimBench --benchmark_filter=ThreadScaling \
 *       --benchmark_out=thread_scaling.json --benchmark_out_format=json
 */

#include <benchmark/benchmark.h>
#include <spdlog/spdlog.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <memory>
#include <vector>

#include "io/StatsWriter.h"
#include "io/checkpointWriter/YAMLWriterCP.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Periodic.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/generators/CuboidGenerator.h"
#include "physics/pairwiseforces/LennardJonesForce.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {

/**
 * @brief Creates domain with periodic boundaries (realistic for large fluid simulations).
 */
Domain createFluidDomain(R3 domain_size, double cutoff) {
    std::array<std::unique_ptr<Boundary>, 6> boundaries;
    boundaries[0] = std::make_unique<Periodic>(BoundaryLocation::LEFT, domain_size, cutoff, 2);
    boundaries[1] = std::make_unique<Periodic>(BoundaryLocation::RIGHT, domain_size, cutoff, 2);
    boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, domain_size, false);
    boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::BACK, domain_size, false);
    boundaries[4] = std::make_unique<Outflow>(BoundaryLocation::UPPER, domain_size);
    boundaries[5] = std::make_unique<Outflow>(BoundaryLocation::LOWER, domain_size);
    return {domain_size, std::move(boundaries)};
}

/**
 * @brief Strong scaling benchmark: Fixed problem size, varying thread count.
 *
 * Tests how well the code parallelizes by keeping particle count constant
 * and varying the number of threads. Ideal scaling would show linear speedup.
 *
 * Measures:
 * - Wall time per iteration
 * - Parallel efficiency = (T_1 / (T_n * n)) * 100%
 * - Speedup = T_1 / T_n
 */
static void bmThreadScalingStrong(benchmark::State& state) {
    spdlog::set_level(spdlog::level::off);

#ifdef _OPENMP
    const int num_threads = static_cast<int>(state.range(0));
    omp_set_num_threads(num_threads);
#else
    if (state.range(0) != 1) {
        state.SkipWithError("OpenMP not enabled, skipping multi-threaded benchmark");
        return;
    }
#endif

    // Fixed problem size: ~10k particles, 50 iterations
    const R3 domain_size = {180.0, 90.0, 1.0};
    const double cutoff = 3.0;
    const size_t num_particles_per_dim = 100;  // 100x100x1 = 10k particles
    const size_t num_iterations = 100;

    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0.0;
    settings.end_time = settings.delta_t * num_iterations;
    settings.cutoff = cutoff;
    settings.dimensions = 2;
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.domain = createFluidDomain(domain_size, cutoff);

    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    auto writer = std::make_unique<XYZWriter>();
    auto cp_writer = std::make_unique<YAMLWriterCP>();
    auto stat_writer = std::make_unique<StatsWriter>();

    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        LinkedCellContainer container(domain_size, cutoff);
        CuboidGenerator generator({10.0, 10.0, 0.0}, {0.0, 0.0, 0.0},
                                  {num_particles_per_dim, num_particles_per_dim, 1U}, {}, 1.0, 1.1225, 0.1, 1.0, 1.0,
                                  0.01);
        ContainerRef particles(container);
        generator.generateParticles(particles, true, false);

        Simulation<LinkedCellContainer> simulation(container, pairwise_forces, single_forces, settings, *writer,
                                                   *cp_writer, *stat_writer);
        state.ResumeTiming();

        simulation.run();

        state.PauseTiming();
        const size_t total_updates = container.size() * num_iterations;
        state.counters["Threads"] = static_cast<double>(num_threads);
        state.counters["Particles"] = static_cast<double>(container.size());
        state.counters["Iterations"] = static_cast<double>(num_iterations);
        state.counters["UpdatesPerSec"] =
            benchmark::Counter(static_cast<double>(total_updates), benchmark::Counter::kIsRate);
        state.ResumeTiming();
    }
}

/**
 * @brief Weak scaling benchmark: Problem size grows with thread count.
 *
 * Tests scalability by increasing particle count proportionally with threads.
 * Ideal weak scaling would show constant time per thread.
 *
 * Each thread gets ~1250 particles to work with.
 */
static void bmThreadScalingWeak(benchmark::State& state) {
    spdlog::set_level(spdlog::level::off);

#ifdef _OPENMP
    const int num_threads = static_cast<int>(state.range(0));
    omp_set_num_threads(num_threads);
#else
    if (state.range(0) != 1) {
        state.SkipWithError("OpenMP not enabled, skipping multi-threaded benchmark");
        return;
    }
#endif

    // Problem size scales with threads: ~1250 particles per thread
    const size_t particles_per_thread = 1250;
    const size_t total_particles = particles_per_thread * static_cast<size_t>(num_threads);
    const auto num_particles_per_dim = static_cast<size_t>(std::sqrt(total_particles));
    const size_t num_iterations = 100;

    // Domain size scales to maintain constant density
    const double base_domain = 50.0;
    const double scale_factor = std::sqrt(static_cast<double>(num_threads));
    const R3 domain_size = {base_domain * scale_factor, base_domain * scale_factor, 1.0};
    const double cutoff = 3.0;

    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0.0;
    settings.end_time = settings.delta_t * num_iterations;
    settings.cutoff = cutoff;
    settings.dimensions = 2;
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.domain = createFluidDomain(domain_size, cutoff);

    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    auto writer = std::make_unique<XYZWriter>();
    auto cp_writer = std::make_unique<YAMLWriterCP>();
    auto stat_writer = std::make_unique<StatsWriter>();

    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        LinkedCellContainer container(domain_size, cutoff);
        CuboidGenerator generator({5.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, {num_particles_per_dim, num_particles_per_dim, 1U},
                                  {}, 1.0, 1.1225, 0.1, 1.0, 1.0, 0.01);
        ContainerRef particles(container);
        generator.generateParticles(particles, true);

        Simulation<LinkedCellContainer> simulation(container, pairwise_forces, single_forces, settings, *writer,
                                                   *cp_writer, *stat_writer);
        state.ResumeTiming();

        simulation.run();

        state.PauseTiming();
        const size_t total_updates = container.size() * num_iterations;
        state.counters["Threads"] = static_cast<double>(num_threads);
        state.counters["Particles"] = static_cast<double>(container.size());
        state.counters["ParticlesPerThread"] = static_cast<double>(container.size()) / num_threads;
        state.counters["Iterations"] = static_cast<double>(num_iterations);
        state.counters["UpdatesPerSec"] =
            benchmark::Counter(static_cast<double>(total_updates), benchmark::Counter::kIsRate);
        state.ResumeTiming();
    }
}

/**
 * @brief Benchmark force calculation phase only (most parallel part).
 */
static void bmThreadScalingForceOnly(benchmark::State& state) {
    spdlog::set_level(spdlog::level::off);

#ifdef _OPENMP
    const int num_threads = static_cast<int>(state.range(0));
    omp_set_num_threads(num_threads);
#else
    if (state.range(0) != 1) {
        state.SkipWithError("OpenMP not enabled, skipping multi-threaded benchmark");
        return;
    }
#endif

    const R3 domain_size = {180.0, 90.0, 1.0};
    const double cutoff = 3.0;
    const size_t num_particles_per_dim = 100;

    LinkedCellContainer container(domain_size, cutoff);
    CuboidGenerator generator({10.0, 10.0, 0.0}, {0.0, 0.0, 0.0}, {num_particles_per_dim, num_particles_per_dim, 1U},
                              {}, 1.0, 1.1225, 0.1, 1.0, 1.0, 0.01);
    ContainerRef particles(container);
    generator.generateParticles(particles, true);

    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0.0;
    settings.end_time = 0.0005;
    settings.cutoff = cutoff;
    settings.dimensions = 2;
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.domain = createFluidDomain(domain_size, cutoff);

    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    auto writer = std::make_unique<XYZWriter>();
    auto cp_writer = std::make_unique<YAMLWriterCP>();
    auto stat_writer = std::make_unique<StatsWriter>();

    Simulation<LinkedCellContainer> simulation(container, pairwise_forces, single_forces, settings, *writer, *cp_writer,
                                               *stat_writer);

    // Warm up and prepare simulation state
    container.prepareForParallelIteration();

    for ([[maybe_unused]] auto _ : state) {
        // Reset forces
        for (auto& p : container) {
            p.getF() = R3{};
        }

        // Benchmark only force calculation
        benchmark::DoNotOptimize(container);
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            simulation.calculateF(0);
        }
        benchmark::ClobberMemory();
    }

    state.counters["Threads"] = static_cast<double>(num_threads);
    state.counters["Particles"] = static_cast<double>(container.size());
}

// Register benchmarks with thread counts: 1, 2, 4, 8, 16
BENCHMARK(bmThreadScalingStrong)
    ->Name("ThreadScaling/Strong/10kParticles")
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->ReportAggregatesOnly(true);

BENCHMARK(bmThreadScalingWeak)
    ->Name("ThreadScaling/Weak/1250PerThread")
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(5)
    ->ReportAggregatesOnly(true);

BENCHMARK(bmThreadScalingForceOnly)
    ->Name("ThreadScaling/ForceCalc/10kParticles")
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Unit(benchmark::kMicrosecond)
    ->Repetitions(10)
    ->ReportAggregatesOnly(true);

}  // namespace mol_sim
