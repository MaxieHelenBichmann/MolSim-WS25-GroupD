#include "physics/targettedforces/TargetForceSource.h"

#include <gtest/gtest.h>

#include "particles/Particle.h"
#include "testingUtils.h"

namespace mol_sim {

/**
 * @brief Test fixture for TargetForceSource tests.
 */
class TargetForceSourceTest : public testing::Test {
   protected:
    double precision = 1e-9;
};

/**
 * @brief Tests that target force is only applied to type 3 particles.
 */
TEST_F(TargetForceSourceTest, AppliesOnlyToType3) {
    R3 direction = {0.0, 0.0, 1.0};
    double magnitude = 0.8;
    size_t max_iterations = 100;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p_target(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);
    Particle p_normal(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 0);
    Particle p_membrane(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    R3 expected_force = {0.0, 0.0, 0.8};
    R3 zero_force = {0.0, 0.0, 0.0};

    EXPECT_R3_NEAR(expected_force, target_force.applyForce(p_target), precision);
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_normal), precision);
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_membrane), precision);
}

/**
 * @brief Tests that target force is only applied to type 4 particles (membrane targets).
 */
TEST_F(TargetForceSourceTest, AppliesOnlyToType4) {
    R3 direction = {1.0, 0.0, 0.0};
    double magnitude = 1.5;
    size_t max_iterations = 50;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p_membrane_target(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 4);
    Particle p_normal(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 0);
    Particle p_membrane(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    R3 expected_force = {1.5, 0.0, 0.0};
    R3 zero_force = {0.0, 0.0, 0.0};

    EXPECT_R3_NEAR(expected_force, target_force.applyForce(p_membrane_target), precision);
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_normal), precision);
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_membrane), precision);
}

/**
 * @brief Tests that force is applied in the correct direction.
 */
TEST_F(TargetForceSourceTest, AppliesInCorrectDirection) {
    R3 direction_z = {0.0, 0.0, 1.0};
    R3 direction_y = {0.0, 1.0, 0.0};
    R3 direction_x = {1.0, 0.0, 0.0};
    double magnitude = 2.0;
    size_t max_iterations = 10;

    TargetForceSource force_z(direction_z, magnitude, max_iterations);
    TargetForceSource force_y(direction_y, magnitude, max_iterations);
    TargetForceSource force_x(direction_x, magnitude, max_iterations);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);

    R3 expected_z = {0.0, 0.0, 2.0};
    R3 expected_y = {0.0, 2.0, 0.0};
    R3 expected_x = {2.0, 0.0, 0.0};

    EXPECT_R3_NEAR(expected_z, force_z.applyForce(p), precision);
    EXPECT_R3_NEAR(expected_y, force_y.applyForce(p), precision);
    EXPECT_R3_NEAR(expected_x, force_x.applyForce(p), precision);
}

/**
 * @brief Tests that force has correct magnitude.
 */
TEST_F(TargetForceSourceTest, CorrectMagnitude) {
    R3 direction = {0.0, 1.0, 0.0};
    double magnitude = 5.5;
    size_t max_iterations = 100;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);

    R3 result = target_force.applyForce(p);
    R3 expected = {0.0, 5.5, 0.0};

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests that force stops after max_iterations.
 */
TEST_F(TargetForceSourceTest, StopsAfterMaxIterations) {
    R3 direction = {0.0, 0.0, 1.0};
    double magnitude = 1.0;
    size_t max_iterations = 3;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);

    R3 expected_force = {0.0, 0.0, 1.0};
    R3 zero_force = {0.0, 0.0, 0.0};

    // First 3 iterations should apply force
    EXPECT_R3_NEAR(expected_force, target_force.applyForce(p), precision);
    EXPECT_R3_NEAR(expected_force, target_force.applyForce(p), precision);
    EXPECT_R3_NEAR(expected_force, target_force.applyForce(p), precision);

    // After max_iterations, should return zero
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p), precision);
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p), precision);
}

/**
 * @brief Tests that iteration counter increments only for target particles.
 */
TEST_F(TargetForceSourceTest, IterationCounterOnlyIncrementsForTargets) {
    R3 direction = {0.0, 0.0, 1.0};
    double magnitude = 1.0;
    size_t max_iterations = 2;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p_target(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);
    Particle p_normal(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 0);

    R3 expected_force = {0.0, 0.0, 1.0};
    R3 zero_force = {0.0, 0.0, 0.0};

    // Apply to normal particle (should not increment counter)
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_normal), precision);
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_normal), precision);
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_normal), precision);

    // Should still have 2 iterations left for target particles
    EXPECT_R3_NEAR(expected_force, target_force.applyForce(p_target), precision);
    EXPECT_R3_NEAR(expected_force, target_force.applyForce(p_target), precision);

    // Now should be exhausted
    EXPECT_R3_NEAR(zero_force, target_force.applyForce(p_target), precision);
}

/**
 * @brief Tests that force works with arbitrary direction vector (not unit vector).
 */
TEST_F(TargetForceSourceTest, ArbitraryDirectionVector) {
    R3 direction = {1.0, 2.0, 3.0};  // Not normalized
    double magnitude = 2.0;
    size_t max_iterations = 10;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 4);

    R3 result = target_force.applyForce(p);
    R3 expected = {2.0, 4.0, 6.0};  // magnitude * direction

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests that force is zero when magnitude is zero.
 */
TEST_F(TargetForceSourceTest, ZeroMagnitude) {
    R3 direction = {0.0, 0.0, 1.0};
    double magnitude = 0.0;
    size_t max_iterations = 100;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);

    R3 result = target_force.applyForce(p);
    R3 expected = {0.0, 0.0, 0.0};

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests getters return correct values.
 */
TEST_F(TargetForceSourceTest, GettersReturnCorrectValues) {
    R3 direction = {1.0, 2.0, 3.0};
    double magnitude = 4.5;
    size_t max_iterations = 150;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    EXPECT_R3_NEAR(direction, target_force.getDirection(), precision);
    EXPECT_NEAR(magnitude, target_force.getMagnitude(), precision);
    EXPECT_EQ(max_iterations, target_force.getMaxIterations());
}

/**
 * @brief Tests that both type 3 and type 4 particles receive the force.
 */
TEST_F(TargetForceSourceTest, AppliesToBothType3And4) {
    R3 direction = {0.0, 1.0, 0.0};
    double magnitude = 0.8;
    size_t max_iterations = 100;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p_type3(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);
    Particle p_type4(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 4);

    R3 expected = {0.0, 0.8, 0.0};

    EXPECT_R3_NEAR(expected, target_force.applyForce(p_type3), precision);
    EXPECT_R3_NEAR(expected, target_force.applyForce(p_type4), precision);
}

/**
 * @brief Tests max_iterations = 0 means no force is applied.
 */
TEST_F(TargetForceSourceTest, ZeroMaxIterations) {
    R3 direction = {0.0, 0.0, 1.0};
    double magnitude = 5.0;
    size_t max_iterations = 0;
    TargetForceSource target_force(direction, magnitude, max_iterations);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 3);

    R3 result = target_force.applyForce(p);
    R3 expected = {0.0, 0.0, 0.0};

    EXPECT_R3_NEAR(expected, result, precision);
}

}  // namespace mol_sim
