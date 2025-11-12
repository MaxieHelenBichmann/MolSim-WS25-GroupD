#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H
#include <random>

#include "particles/Particle.h"

namespace mol_sim {
/**
 * @brief Creates a random R3 Vector
 *
 * @return R3
 */
inline R3 randomR3() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> dist{};

    return {dist(gen), dist(gen), dist(gen)};
}

/**
 * @brief Creates a random particle
 *
 * @return Particle
 */

inline Particle randomParticle() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> dist{};

    return {randomR3(), randomR3(), dist(gen), dist(gen), dist(gen)};
}

}  // namespace mol_sim
#endif
