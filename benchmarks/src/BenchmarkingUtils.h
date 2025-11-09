#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H
#include <random>

#include "particles/Particle.h"
#include "spdlog/spdlog.h"

namespace mol_sim {
inline R3 randomR3() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> dist{};

    return {dist(gen), dist(gen), dist(gen)};
}

inline Particle randomParticle() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> dist{};

    return {randomR3(), randomR3(), dist(gen), dist(gen), dist(gen)};
}

inline void setUpLogging() { spdlog::set_level(spdlog::level::warn); }
}  // namespace mol_sim
#endif
