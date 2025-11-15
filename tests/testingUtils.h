#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include <gtest/gtest.h>

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
inline ::testing::AssertionResult equalityR3(const R3& expected, const R3& actual, double tolerance = 1e-9) {
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

#define EXPECT_R3_EQUAL(expected, actual) EXPECT_TRUE(AreR3sEqual(expected, actual))

#define EXPECT_R3_EQUAL_WITH_TOLERANCE(expected, actual, tolerance) \
    EXPECT_TRUE(AreR3sEqual(expected, actual, tolerance))
}  // namespace mol_sim

#endif
