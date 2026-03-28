#include <gtest/gtest.h>

#include <numbers>

#include "particles/Particle.h"
#include "physics/pairwiseforces/LennardJonesForce.h"
#include "testingUtils.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the Lennard-Jones force
 * Base data/setup:
 * p1:
 *  position: (0,0,0)
 *  velocity: (0,0,0)
 *  (mass: 1)
 * p2:
 *  position: individually set in tests
 *  velocity: (0,0,0)
 *  (mass: 1)
 * big_epsilon:
 *  p1 has a non-shallow potential well (epsilon = 5)
 * small_epsilon
 *  p1 has a shallow potential well (epsilon = 1)
 * big_sigma:
 *  p2 replused already at large distance to p1 (sigma = 5)
 * small_sigma
 *  p2 can get relatively close to p1 (sigma = 1)
 */
class LennardJonesForceTest : public testing::Test {
   protected:
    R3 p1_x = {0.0, 0.0, 0.0};
    R3 p2_x = {0.0, 0.0, 0.0};
    R3 p1_v = {0.0, 0.0, 0.0};
    R3 p2_v = {0.0, 0.0, 0.0};

    Particle p1;
    Particle p2;
    double big_epsilon = 5;
    double small_epsilon = 1;
    double big_sigma = 5;
    double small_sigma = 1;
    double precision = 5e-8;

    LennardJonesForceTest() : p1{p1_x, p1_v, 1.0, 5., 1.}, p2{p2_x, p2_v, 1.0, 5., 1.} {}

    void SetUp() override {
        p1.getEpsilon() = big_epsilon;
        p2.getEpsilon() = big_epsilon;
        p1.getSigma() = small_sigma;
        p2.getSigma() = small_sigma;
    }
};

/**
 * @brief Tests the Lennard Jones Force between two far away particles that are aligned on the x axis.
 * Compares against hand computed values.
 *
 * */
TEST_F(LennardJonesForceTest, FarAway_X_Aligned) {
    R3 expected_force = {0.007320642471, 0., 0.};
    p2.getX() = {4.0, 0.0, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
    p2.getX() = {-4.0, 0.0, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}
/**
 * @brief Tests the Lennard Jones Force between two far away particles that are aligned on the y axis.
 * Compares against hand computed values.
 *
 * */
TEST_F(LennardJonesForceTest, FarAway_Y_Aligned) {
    R3 expected_force = {0., 0.007320642471, 0.};
    p2.getX() = {0.0, 4.0, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);

    EXPECT_R3_NEAR(expected_force, result, precision);

    p2.getX() = {0.0, -4.0, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);

    EXPECT_R3_NEAR(-1 * expected_force, result, precision);
}
/**
 * @brief Tests the Lennard Jones Force between two far away particles that are not aligned.
 * Compares against hand computed values.
 *
 * */
TEST_F(LennardJonesForceTest, FarAway_Diagonal) {
    R3 expected_force = {0.0004577357322, 0.0004577357322, 0.};
    p2.getX() = {4.0, 4.0, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);

    p2.getX() = {-4.0, -4.0, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(-1 * expected_force, result, precision);
}

/**
 * @brief Tests the computed force against hand computed force when particle is at bottom of potential
 * well. Force should be 0 there.
 *
 * */
class LennardJonesForceAtBottomTest : public LennardJonesForceTest {
   protected:
    void SetUp() override { LennardJonesForceTest::SetUp(); }
};

/**
 * @brief Testing the afromentioned case at radius 2 when the two particles are cardinal to one another.
 *
 */

TEST_F(LennardJonesForceAtBottomTest, PotentialWellAtRadius2_Cardinal) {
    double sigma = 1.78179743599999995673499597615;  // 2^{5/6}
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    const double expected_force = 0.0;

    p2.getX() = {2.0, 0.0, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_NEAR(result[0], expected_force, precision);
    EXPECT_EQ(result[1], 0.);
    EXPECT_EQ(result[2], 0.);

    p2.getX() = {0.0, -2.0, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_EQ(result[0], 0);
    EXPECT_NEAR(result[1], expected_force, precision);
    EXPECT_EQ(result[2], 0.);
}

/**
 * @brief Testing the afromentioned case at radius 2 when the two particles are diagonal to one another.
 *
 */

TEST_F(LennardJonesForceAtBottomTest, PotentialWellAtRadius2_Diagonal) {
    double sigma = 1.78179743599999995673499597615;  // 2^{5/6}
    double bottom_diagonal = std::numbers::sqrt2;
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    R3 expected_force = {0., 0., 0.};

    p2.getX() = {bottom_diagonal, bottom_diagonal, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
    p2.getX() = {-bottom_diagonal, -bottom_diagonal, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
}
/**
 * @brief Testing the afromentioned case at radius 1 when the two particles are cardinal to one another.
 *
 */

TEST_F(LennardJonesForceAtBottomTest, PotentialWellAtRadius1_Cardinal) {
    double sigma = 8.90898718099999986641535087983E-1;  // 2^{-1/6}
    p1.getEpsilon() = small_epsilon;
    p2.getEpsilon() = small_epsilon;
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    R3 expected_force = {0., 0., 0.};

    p2.getX() = {1.0, 0.0, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);

    p2.getX() = {0.0, 1.0, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Testing the afromentioned case at radius 1 when the two particles are diagonal to one another.
 *
 */

TEST_F(LennardJonesForceAtBottomTest, PotentialWellAtRadius1_Diagonal) {
    double sigma = 8.90898718099999986641535087983E-1;  // 2^{-1/6}
    double bottom_diagonal = 1.0 / std::numbers::sqrt2;
    p1.getEpsilon() = small_epsilon;
    p2.getEpsilon() = small_epsilon;
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    R3 expected_force = {0., 0., 0.};

    p2.getX() = {bottom_diagonal, bottom_diagonal, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
    p2.getX() = {bottom_diagonal, -bottom_diagonal, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Tests compute force against hand computed force for particle thats exactly at
 * sigma of the other particle's potential well (i.e. at the
 * zero crossing of the potential). We should experience a strong
 * repulsion at this point.
 *
 * */
class LennardJonesForceAtSigmaTest : public LennardJonesForceTest {
   protected:
    void SetUp() override { LennardJonesForceTest::SetUp(); }
};

/**
 * @brief Testing the afromentioned case when the two particles are cardinal to one another.
 *
 */
TEST_F(LennardJonesForceAtSigmaTest, Cardinal) {
    R3 expected_force = {-120.0, 0., 0.};
    p2.getX() = {small_sigma, 0.0, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);

    p2.getX() = {-small_sigma, 0.0, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(-1 * expected_force, result, precision);

    p2.getX() = {0.0, small_sigma, 0.0};
    expected_force = {0., -120.0, 0.};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
}
/**
 * @brief Testing the afromentioned case when the two particles are diagonal to one another.
 *
 */
TEST_F(LennardJonesForceAtSigmaTest, Diagonal) {
    R3 expected_force = {-60.0 * std::numbers::sqrt2, -60.0 * std::numbers::sqrt2, 0.};
    double sigma_diagonal = 1.0 / std::numbers::sqrt2;

    p2.getX() = {sigma_diagonal, sigma_diagonal, 0.0};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);

    expected_force = {-60.0 * std::numbers::sqrt2, 60.0 * std::numbers::sqrt2, 0.};
    p2.getX() = {sigma_diagonal, -sigma_diagonal, 0.0};
    result = LennardJonesForce().applyForce(p1, p2);
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Tests that Newton's 3rd law holds for our force calculations.
 * (i.e. force particle 1 = -force particle 2).
 *
 * */
TEST_F(LennardJonesForceTest, Newton3) {
    p1.getEpsilon() = big_epsilon;
    p2.getEpsilon() = big_epsilon;
    p1.getSigma() = small_sigma;
    p2.getSigma() = small_sigma;
    p2.getX() = {4.0, 4.0, 0.0};
    R3 result1 = LennardJonesForce().applyForce(p1, p2);
    R3 result2 = LennardJonesForce().applyForce(p2, p1);
    EXPECT_R3_NEAR(-1 * result2, result1, precision);
}

/**
 * @brief Tests that the force is 0 when the particles are at the same position.
 *
 */
TEST_F(LennardJonesForceTest, TestForceCalcZeroDistance) {
    p2 = {p1_x, p2_v, 1.0, 5., 1.};
    R3 result = LennardJonesForce().applyForce(p1, p2);
    R3 expected_force = {0., 0., 0.};
    EXPECT_R3_NEAR(expected_force, result, precision);
}
/**
 * @brief Tests that the force between particles with different epsilons and sigmas gets computed correctly.
 *
 */
TEST_F(LennardJonesForceTest, MixingRules) {
    p1.getSigma() = 1.0;
    p1.getEpsilon() = 5.0;
    p2.getSigma() = 2.0;
    p2.getEpsilon() = 10.0;

    p2.getX() = {3.0, 4.0, 0.0};

    const R3 expected_force = {1.48242203050000000141839251455E-2, 1.97656040600000006135505969951E-2, 0.};
    R3 result = LennardJonesForce().applyForce(p1, p2);

    EXPECT_R3_NEAR(expected_force, result, precision);
}

}  // namespace mol_sim
