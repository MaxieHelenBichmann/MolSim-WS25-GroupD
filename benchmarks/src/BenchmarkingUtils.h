#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H
#include <array>
#include <random>

#include "particles/Particle.h"

namespace mol_sim {
inline std::mt19937& getGenerator() {
    static thread_local std::random_device rd{};
    static thread_local std::mt19937 gen{rd()};
    return gen;
}

inline double randomBounded(double lower, double upper) {
    std::normal_distribution<double> dist{};
    double res;
    do {
        res = dist(getGenerator());
    } while (res < lower || res > upper);
    return res;
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

}  // namespace mol_sim
#endif
