#include "physics/singleforces/GravForce.h"

#include <gtest/gtest.h>

#include "particles/Particle.h"
#include "testingUtils.h"

namespace mol_sim {

/**
 * @brief Test fixture for GravForce tests.
 */
class GravForceTest : public testing::Test {
   protected:
    double precision = 1e-9;
};

/**
 * @brief Tests gravitational force in negative z direction (downward).
 */
TEST_F(GravForceTest, GravityDownwardZ) {
    R3 g_grav = {0.0, 0.0, -9.81};
    GravForce grav_force(g_grav);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 2.0, 1.0, 1.0, 0);

    R3 result = grav_force.applyForce(p);
    R3 expected = {0.0, 0.0, -19.62};

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests gravitational force in negative y direction (downward).
 */
TEST_F(GravForceTest, GravityDownwardY) {
    R3 g_grav = {0.0, -9.81, 0.0};
    GravForce grav_force(g_grav);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 3.0, 1.0, 1.0, 0);

    R3 result = grav_force.applyForce(p);
    R3 expected = {0.0, -29.43, 0.0};

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests gravitational force in positive x direction.
 */
TEST_F(GravForceTest, GravityPositiveX) {
    R3 g_grav = {5.0, 0.0, 0.0};
    GravForce grav_force(g_grav);

    Particle p(R3{1.0, 1.0, 1.0}, R3{0.0, 0.0, 0.0}, 1.5, 1.0, 1.0, 0);

    R3 result = grav_force.applyForce(p);
    R3 expected = {7.5, 0.0, 0.0};

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests gravitational force with arbitrary direction.
 */
TEST_F(GravForceTest, GravityArbitraryDirection) {
    R3 g_grav = {1.0, -2.0, 3.0};
    GravForce grav_force(g_grav);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 4.0, 1.0, 1.0, 0);

    R3 result = grav_force.applyForce(p);
    R3 expected = {4.0, -8.0, 12.0};

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests that force is zero when gravitational acceleration is zero.
 */
TEST_F(GravForceTest, ZeroGravity) {
    R3 g_grav = {0.0, 0.0, 0.0};
    GravForce grav_force(g_grav);

    Particle p(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 10.0, 1.0, 1.0, 0);

    R3 result = grav_force.applyForce(p);
    R3 expected = {0.0, 0.0, 0.0};

    EXPECT_R3_NEAR(expected, result, precision);
}

/**
 * @brief Tests that force scales correctly with mass.
 */
TEST_F(GravForceTest, ForceMassScaling) {
    R3 g_grav = {0.0, 0.0, -0.001};
    GravForce grav_force(g_grav);

    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 0);
    Particle p2(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 100.0, 1.0, 1.0, 0);

    R3 result1 = grav_force.applyForce(p1);
    R3 result2 = grav_force.applyForce(p2);

    R3 expected1 = {0.0, 0.0, -0.001};
    R3 expected2 = {0.0, 0.0, -0.1};

    EXPECT_R3_NEAR(expected1, result1, precision);
    EXPECT_R3_NEAR(expected2, result2, precision);
}

/**
 * @brief Tests that gravity works regardless of particle type.
 */
TEST_F(GravForceTest, WorksForAllParticleTypes) {
    R3 g_grav = {0.0, -1.0, 0.0};
    GravForce grav_force(g_grav);

    Particle p_normal(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 2.0, 1.0, 1.0, 0);
    Particle p_membrane(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 2.0, 1.0, 1.0, 2);
    Particle p_target(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 2.0, 1.0, 1.0, 3);

    R3 expected = {0.0, -2.0, 0.0};

    EXPECT_R3_NEAR(expected, grav_force.applyForce(p_normal), precision);
    EXPECT_R3_NEAR(expected, grav_force.applyForce(p_membrane), precision);
    EXPECT_R3_NEAR(expected, grav_force.applyForce(p_target), precision);
}

/**
 * @brief Tests getType() returns GRAV.
 */
TEST_F(GravForceTest, ReturnsCorrectType) {
    R3 g_grav = {0.0, -9.81, 0.0};
    GravForce grav_force(g_grav);

    EXPECT_EQ(SingleForce::GRAV, grav_force.getType());
}

}  // namespace mol_sim
