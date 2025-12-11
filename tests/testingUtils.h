#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include <gtest/gtest.h>
#include <particles/Particle.h>

#include <cmath>

#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Compares two R3 objects for equality within a given tolerance.
 *
 * @param expected The expected R3 value.
 * @param actual The actual R3 value.
 * @param tolerance The tolerance for floating-point comparison.
 * @return ::testing::AssertionSuccess() if they are equal, otherwise ::testing::AssertionFailure().
 */
inline ::testing::AssertionResult equalityR3(const R3& actual, const R3& expected, double tolerance = 1e-9) {
    if (std::abs(expected[0] - actual[0]) < tolerance && std::abs(expected[1] - actual[1]) < tolerance &&
        std::abs(expected[2] - actual[2]) < tolerance) {
        return ::testing::AssertionSuccess();
    }

    return ::testing::AssertionFailure() << "R3 objects are not equal within tolerance " << tolerance << ".\n"
                                         << "  Expected: {" << expected[0] << ", " << expected[1] << ", " << expected[2]
                                         << "}\n"
                                         << "  Actual:   {" << actual[0] << ", " << actual[1] << ", " << actual[2]
                                         << "}";
}
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define EXPECT_R3_EQ(expected, actual) EXPECT_TRUE(equalityR3(expected, actual))
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define EXPECT_R3_NEAR(expected, actual, tolerance) EXPECT_TRUE(equalityR3(expected, actual, tolerance))

inline ::testing::AssertionResult equalityParticle(const mol_sim::Particle& actual, const mol_sim::Particle& expected,
                                                   double tolerance = 1e-9) {
    if (actual.getEpsilon() == expected.getEpsilon() && actual.getM() == expected.getM() &&
        actual.getSigma() == expected.getSigma() && actual.getType() == expected.getType() &&
        equalityR3(actual.getF(), expected.getF(), tolerance) &&
        equalityR3(actual.getOldF(), expected.getOldF(), tolerance) &&
        equalityR3(actual.getV(), expected.getV(), tolerance) &&
        equalityR3(actual.getX(), expected.getX(), tolerance) &&
        equalityR3(actual.getOldX(), expected.getOldX(), tolerance)) {
        return ::testing::AssertionSuccess();
    }

    return ::testing::AssertionFailure() << "R3 objects are not equal within tolerance " << tolerance << ".\n"
                                         << "  Expected: " << expected.toString() << "\n"
                                         << "  Actual: " << actual.toString() << "}";
}
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define EXPECT_PARTICLE_EQ(expected, actual) EXPECT_TRUE(equalityParticle(expected, actual))
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define EXPECT_PARTICLE_NEAR(expected, actual, tolerance) EXPECT_TRUE(equalityParticle(expected, actual, tolerance))
}  // namespace mol_sim
#endif
