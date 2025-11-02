#include <gtest/gtest.h>

#include "particles/Particle.h"
#include "physics/GravitationalForce.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the GravitationalForce
 * Base data/setup:
 * p1:
 *  position: (0,0,0)
 *  velocity: (0,0,0)
 *  mass: 1
 * p2:
 *  position: (1,0,0)
 *  velocity: (0,0,0)
 *  mass: 2
 */
class GravitationalForceTest : public testing::Test {
   protected:
    R3 p1_x = {1.0, 0.5, 2.0};
    R3 p2_x = {1.0, 2.0, 4.0};
    R3 p1_v = {0.0, 0.0, 0.0};
    R3 p2_v = {0.0, 0.0, 0.0};
    double p1_mass = 1.0;
    double p2_mass = 2.0;

    Particle p1;
    Particle p2;

    GravitationalForceTest() : p1{p1_x, p1_v, p1_mass}, p2{p2_x, p2_v, p2_mass} {}
};

/**
 * @brief Tests that the force is computed correctly when both particles are not axis aligned.
 *
 */
TEST_F(GravitationalForceTest, TestForceCalcUnaligned) {
    R3 result = GravitationalForce().calculateForce(p1, p2);
    R3 expected = {0.0, 0.192, 0.256};
    EXPECT_EQ(result, expected);
}
/**
 * @brief Tests that in force calculations Newtons third law holds.
 * Force on particle1 = -Force on particle2
 */
TEST_F(GravitationalForceTest, TestForceSymmetry) {
    R3 result1 = GravitationalForce().calculateForce(p1, p2);
    R3 result2 = GravitationalForce().calculateForce(p2, p1);
    EXPECT_EQ(result1, result2 * -1);
}

/**
 * @brief Tests that the force is 0 when the particles are at the same position.
 *
 */
TEST_F(GravitationalForceTest, TestForceCalcZeroDistance) {
    p2 = {p1_x, p2_v, p2_mass};
    R3 result = GravitationalForce().calculateForce(p1, p2);
    R3 expected = {0.0, 0.0, 0.0};
    EXPECT_EQ(result, expected);
}
/**
 * @brief Tests that the force is 0 when one of the particles has mass 0.
 *
 */
TEST_F(GravitationalForceTest, TestForceCalcZeroMass) {
    p1 = {p1_x, p1_v, 0.0};
    R3 result = GravitationalForce().calculateForce(p1, p2);
    R3 expected = {0.0, 0.0, 0.0};
    EXPECT_EQ(result, expected);
}
/**
 * @brief Tests that the force is computed correctly when both particles are axis aligned.
 *
 */
TEST_F(GravitationalForceTest, TestForceCalcAligned) {
    p1 = {{0.0, 0.0, 0.0}, p1_v, 1.0};
    p2 = {{1.0, 0.0, 0.0}, p2_v, 1.0};
    R3 result = GravitationalForce().calculateForce(p1, p2);
    R3 expected = {1.0, 0.0, 0.0};
    EXPECT_EQ(result, expected);
}

}  // namespace mol_sim
