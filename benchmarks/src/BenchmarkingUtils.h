#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H
#include <array>
#include <memory>
#include <random>

#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Periodic.h"
#include "particles/boundaries/Reflecting.h"
#include "utils/Settings.h"

namespace mol_sim {
inline std::mt19937& getGenerator() {
    static thread_local std::random_device rd{};
    static thread_local std::mt19937 gen{rd()};
    return gen;
}

inline double randomBounded(double lower, double upper) {
    std::uniform_real_distribution<double> dist{lower, upper};
    return dist(getGenerator());
}
/**
 * @brief Creates a random R3 Vector
 *
 * @return R3
 */
inline R3 randomR3() {
    std::normal_distribution<double> dist{};
    auto& gen = getGenerator();
    return {dist(gen), dist(gen), dist(gen)};
}
inline R3 randomR3(std::array<double, 6> bounds) {
    double x = randomBounded(bounds[0], bounds[1]);
    double y = randomBounded(bounds[2], bounds[3]);
    double z = randomBounded(bounds[4], bounds[5]);

    return {x, y, z};
}

/**
 * @brief Creates a random particle
 *
 * @return Particle
 */

inline Particle randomParticle() {
    std::normal_distribution<double> dist{};
    auto& gen = getGenerator();
    return {randomR3(), randomR3(), dist(gen), dist(gen), dist(gen)};
}

/**
 * 2D domain with periodic boundaries on left/right and reflecting boundaries on top/bottom.
 */
inline SettingsParam createContestSettings() {
    SettingsParam settings;
    settings.delta_t = 0.0005;
    settings.start_time = 0.0;
    settings.end_time = 0.5;
    settings.cutoff = 3;
    settings.dimensions = 2;
    settings.base_name = "contest";
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.single_forces = {SingleForce::GRAV};
    settings.thermo = true;
    settings.init_temp = 40;
    settings.target_temp = 40.;
    settings.thermostat_freq = 1000;
    settings.g_grav_vec = {0.0, -12.44, 0.0};
    R3 domain_size = {300., 54., 1.};
    std::array<std::unique_ptr<Boundary>, 6> boundaries;
    boundaries[0] = std::make_unique<Periodic>(BoundaryLocation::LEFT, domain_size, 3, 2);
    boundaries[1] = std::make_unique<Periodic>(BoundaryLocation::RIGHT, domain_size, 3, 2);
    boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, domain_size, false);
    boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::BACK, domain_size, false);
    boundaries[4] = std::make_unique<Outflow>(BoundaryLocation::UPPER, domain_size);
    boundaries[5] = std::make_unique<Outflow>(BoundaryLocation::LOWER, domain_size);

    settings.domain = Domain(domain_size, std::move(boundaries));
    return settings;
}

}  // namespace mol_sim
#endif
