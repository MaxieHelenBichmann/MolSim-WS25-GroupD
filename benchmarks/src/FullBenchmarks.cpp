/**
 * @file FullBenchmarks.cpp
 * @brief Performance benchmarks for large-scale molecular simulations.
 *
 * Measures runtime and molecule-updates-per-second for representative scenarios with ≥10000 particles.
 * Tests 2D simulation with periodic boundaries (left/right) and reflecting boundaries (top/bottom).
 * I/O is disabled to measure pure simulation performance.
 *
 * Run with: ./build/benchmarks/benchmarks --benchmark_filter=FullSimulation
 */
#include <benchmark/benchmark.h>
#include <spdlog/spdlog.h>

#include <array>
#include <memory>

#include "io/checkpointWriter/YAMLWriterCP.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Periodic.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/container/LinkedCellContainer.h"
#include "physics/LennardJonesForce.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {

/**
 * @brief Generates a cuboid of particles directly into any ParticleContainer.
 */
template <ParticleContainer Container>
void generateCuboid(Container& particles, R3 position, R3 velocity, Vector<size_t, 3> num_particles, double mass,
                    double distance, double avg_velo, double epsilon, double sigma) {
    particles.reserve(particles.size() + (num_particles[0] * num_particles[1] * num_particles[2]));
    for (size_t i = 0; i < num_particles[2]; i++) {
        for (size_t j = 0; j < num_particles[1]; j++) {
            for (size_t k = 0; k < num_particles[0]; k++) {
                R3 curr_pos = {position[0] + (static_cast<double>(k) * distance),
                               position[1] + (static_cast<double>(j) * distance),
                               position[2] + (static_cast<double>(i) * distance)};
                R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
                particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma);
            }
        }
    }
}

/**
 * @brief Creates simulation settings for large-scale benchmark.
 * 2D domain with periodic boundaries on left/right and reflecting boundaries on top/bottom.
 */
SettingsParam createBenchmarkSettings(R3 domain_size, double cutoff, double delta_t, double end_time) {
    SettingsParam settings;
    settings.delta_t = delta_t;
    settings.start_time = 0.0;
    settings.end_time = end_time;
    settings.epsilon = 1.0;
    settings.sigma = 1.0;
    settings.cutoff = cutoff;
    settings.dimensions = 2;
    settings.base_name = "benchmark";
    settings.force = Force::LENNARDJONES;
    settings.thermo = true;
    settings.init_temp = 20;
    settings.target_temp = 40.;
    settings.delta_temp = 0.1;
    settings.thermostat_freq = 1000;
    settings.g_grav = -12.44;

    std::array<std::unique_ptr<Boundary>, 6> boundaries;
    boundaries[0] = std::make_unique<Periodic>(BoundaryLocation::LEFT, domain_size, cutoff, 2);
    boundaries[1] = std::make_unique<Periodic>(BoundaryLocation::RIGHT, domain_size, cutoff, 2);
    boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, domain_size, false);
    boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::BACK, domain_size, false);
    boundaries[4] = std::make_unique<Outflow>(BoundaryLocation::UPPER, domain_size);
    boundaries[5] = std::make_unique<Outflow>(BoundaryLocation::LOWER, domain_size);

    settings.domain = Domain(domain_size, std::move(boundaries));
    return settings;
}

/**
 * @brief Benchmarks full simulation loop with LinkedCellContainer (10000+ particles).
 * Setup: 100x100x1 particle grid in 2D domain with periodic left/right, reflecting top/bottom.
 * Measures: Total runtime and calculates molecule-updates-per-second.
 */
static void bmSimulationFull(benchmark::State& state) {
    spdlog::set_level(spdlog::level::warn);

    const R3 domain_size = {150.0, 150.0, 1.0};
    const double cutoff = 3.0;
    const double delta_t = 0.0005;
    const double end_time = 3.0;

    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    auto cp_writer = std::make_unique<YAMLWriterCP>();

    for ([[maybe_unused]] auto _ : state) {
        state.PauseTiming();
        LinkedCellContainer container(domain_size, cutoff);
        SettingsParam settings = createBenchmarkSettings(domain_size, cutoff, delta_t, end_time);

        generateCuboid(container, {5.0, .5, 0.0}, {0.0, -25.0, 0.0}, {100U, 100U, 1U}, 1.0, 1.2, 0.1, 1.0, 1.0);

        const size_t num_particles = container.size();
        const auto num_iterations = static_cast<size_t>((end_time - settings.start_time) / delta_t);
        state.counters["Particles"] = static_cast<double>(num_particles);
        state.counters["Iterations"] = static_cast<double>(num_iterations);

        Simulation<LinkedCellContainer> simulation(container, *force_source, settings, *writer, *cp_writer);
        state.ResumeTiming();

        simulation.run();

        state.PauseTiming();
        const size_t total_updates = num_particles * num_iterations;
        state.counters["MoleculeUpdates"] = static_cast<double>(total_updates);
        state.counters["UpdatesPerSec"] =
            benchmark::Counter(static_cast<double>(total_updates), benchmark::Counter::kIsRate);
        state.ResumeTiming();
    }
}

BENCHMARK(bmSimulationFull)->Name("Simulation/Full/LinkedCell")->Unit(benchmark::kMillisecond)->Repetitions(5);

}  // namespace mol_sim
