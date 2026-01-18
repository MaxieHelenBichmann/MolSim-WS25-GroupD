#include <gtest/gtest.h>

#include <numbers>

#include "particles/Particle.h"
#include "physics/SmoothLennardJonesForce.h"
#include "testingUtils.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the smooth Lennard-Jones force
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
class SmoothLennardJonesForceTest : public testing::Test {
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

    double cutoff_radius = 30.0;
    double smoothing_radius = 10.0;

    SmoothLennardJonesForce force;

    SmoothLennardJonesForceTest() : p1{p1_x, p1_v, 1.0, 5., 1.}, p2{p2_x, p2_v, 1.0, 5., 1.} {}

    void SetUp() override {
        p1.getEpsilon() = big_epsilon;
        p2.getEpsilon() = big_epsilon;
        p1.getSigma() = small_sigma;
        p2.getSigma() = small_sigma;
        SmoothLennardJonesForce().initForce(cutoff_radius, smoothing_radius);
    }
};

// ----------------------------------------------------------------------------------------------------
/**
 * @brief Tests the Smooth Lennard Jones Force between two far away particles that are not aligned.
 *
 *  Should be zero due to explicit cutoff.
 */
class SmoothLennardJonesFarAwayTest : public SmoothLennardJonesForceTest {
   protected:
    void SetUp() override { SmoothLennardJonesForceTest::SetUp(); }
};

/**
 * @brief Tests the Smooth Lennard Jones Force between two far away particles that are aligned on the x axis.
 *
 * */
TEST_F(SmoothLennardJonesFarAwayTest, FarAway_Cardinal) {
    R3 expected_force = {0., 0., 0.};
    p2.getX() = {100.0, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
    p2.getX() = {-100.0, 0.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}
/**
 * @brief Tests the Smooth Lennard Jones Force between two far away particles that are not aligned.
 *
 * */
TEST_F(SmoothLennardJonesFarAwayTest, FarAway_Diagonal) {
    R3 expected_force = {0., 0., 0.};
    p2.getX() = {100.0, 100.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-100.0, -100.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}

// ----------------------------------------------------------------------------------------------------
/**
 * @brief Tests the Smooth Lennard Jones Force between two close particles.
 *
 *  Should be exactly equal to the normal Lennard Jones Force.
 */
class SmoothLennardJonesCloseTest : public SmoothLennardJonesForceTest {
   protected:
    void SetUp() override { SmoothLennardJonesForceTest::SetUp(); }
};

/**
 * @brief Tests the Smooth Lennard Jones Force between close particles that are not aligned.
 * Compares against hand computed values.
 *
 * */
TEST_F(SmoothLennardJonesCloseTest, Close_LJ) {
    R3 expected_force = {0.0004577357322, 0.0004577357322, 0.};
    p2.getX() = {4.0, 4.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-4.0, -4.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}

/**
 * @brief Testing the case when the two close particles are in the bottom well (no force), when they are cardinal to one
 * another.
 */
TEST_F(SmoothLennardJonesCloseTest, Close_PotentialWellAtRadius2_Cardinal) {
    double sigma = 1.78179743599999995673499597615;  // 2^{5/6}
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    const double expected_force = 0.0;

    p2.getX() = {2.0, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_NEAR(result[0], expected_force, precision);
    EXPECT_EQ(result[1], 0.);
    EXPECT_EQ(result[2], 0.);

    p2.getX() = {0.0, -2.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_EQ(result[0], 0);
    EXPECT_NEAR(result[1], expected_force, precision);
    EXPECT_EQ(result[2], 0.);
}

/**
 * @brief Testing the case when the two close particles are in the bottom well (no force), when they are diagonal to one
 * another.
 */
TEST_F(SmoothLennardJonesCloseTest, Close_PotentialWellAtRadius2_Diagonal) {
    double sigma = 1.78179743599999995673499597615;  // 2^{5/6}
    double bottom_diagonal = std::numbers::sqrt2;
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    R3 expected_force = {0., 0., 0.};

    p2.getX() = {bottom_diagonal, bottom_diagonal, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
    p2.getX() = {-bottom_diagonal, -bottom_diagonal, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
}
/**
 * @brief Testing the case when the two close particles are in the bottom well (no force), when they are cardinal to one
 * another.
 */
TEST_F(SmoothLennardJonesCloseTest, Close_PotentialWellAtRadius1_Cardinal) {
    double sigma = 8.90898718099999986641535087983E-1;  // 2^{-1/6}
    p1.getEpsilon() = small_epsilon;
    p2.getEpsilon() = small_epsilon;
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    R3 expected_force = {0., 0., 0.};

    p2.getX() = {1.0, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {0.0, 1.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
}

/**
 * @brief Testing the case when the two close particles are in the bottom well (no force), when they are diagonal to one
 * another.
 */
TEST_F(SmoothLennardJonesCloseTest, Close_PotentialWellAtRadius1_Diagonal) {
    double sigma = 8.90898718099999986641535087983E-1;  // 2^{-1/6}
    double bottom_diagonal = 1.0 / std::numbers::sqrt2;
    p1.getEpsilon() = small_epsilon;
    p2.getEpsilon() = small_epsilon;
    p1.getSigma() = sigma;
    p2.getSigma() = sigma;
    R3 expected_force = {0., 0., 0.};

    p2.getX() = {bottom_diagonal, bottom_diagonal, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
    p2.getX() = {bottom_diagonal, -bottom_diagonal, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
}

/**
 * @brief Testing the case when a particle is exactly at the sigma of the other particle's potential well (strong
 * repulsion), and the two particles are cardinal to one another.
 */
TEST_F(SmoothLennardJonesCloseTest, Close_Cardinal) {
    R3 expected_force = {-120.0, 0., 0.};
    p2.getX() = {small_sigma, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-small_sigma, 0.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);

    p2.getX() = {0.0, small_sigma, 0.0};
    expected_force = {0., -120.0, 0.};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
}
/**
 * @brief Testing the case when a particle is exactly at the sigma of the other particle's potential well (strong
 * repulsion), and the two particles are diagonal to one another.
 */
TEST_F(SmoothLennardJonesCloseTest, Close_Diagonal) {
    R3 expected_force = {-60.0 * std::numbers::sqrt2, -60.0 * std::numbers::sqrt2, 0.};
    double sigma_diagonal = 1.0 / std::numbers::sqrt2;

    p2.getX() = {sigma_diagonal, sigma_diagonal, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    expected_force = {-60.0 * std::numbers::sqrt2, 60.0 * std::numbers::sqrt2, 0.};
    p2.getX() = {sigma_diagonal, -sigma_diagonal, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);
}
// ----------------------------------------------------------------------------------------------------
/**
 * @brief Tests the Smooth Lennard Jones Force between two particles at a medium distance, meaning between the smoothing
 * and cutoff radius.
 */
class SmoothLennardJonesMediumDistanceTest : public SmoothLennardJonesForceTest {
   protected:
    void SetUp() override { SmoothLennardJonesForceTest::SetUp(); }
};

/**
 * @brief Test the Smooth Lennard Jones Force between medium distance particles that are aligned on the x axis.
 * Compares against hand computed values.
 */
TEST_F(SmoothLennardJonesMediumDistanceTest, MediumDistance_Cardinal) {
    SmoothLennardJonesForce force;
    force.initForce(cutoff_radius, 3);

    R3 expected_force = {0.018655817750600174, 0., 0.};
    p2.getX() = {3.5, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-3.5, 0.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}

/**
 * @brief Test the Smooth Lennard Jones Force between medium distance particles that are diagonal.
 */
TEST_F(SmoothLennardJonesMediumDistanceTest, MediumDistance_Diagonal) {
    SmoothLennardJonesForce force;
    force.initForce(cutoff_radius, 3);

    R3 expected_force = {0.00026160748585069596, 0.0007848224575520879, 0.0010464299434027838};
    p2.getX() = {1.0, 3.0, 4.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-1.0, -3.0, -4.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}

/**
 * @brief Test the Smooth Lennard Jones Force between medium distance particles when all smoothing is applied
 * (i.e. smoothing radius = 0).
 */
TEST_F(SmoothLennardJonesMediumDistanceTest, MediumDistance_AllSmoothing) {
    SmoothLennardJonesForce force;
    force.initForce(cutoff_radius, 0);

    R3 expected_force = {8.327560625865901, 0.0, 0.0};
    p2.getX() = {1.4, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-1.4, 0.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}
/**
 * @brief Test the Smooth Lennard Jones Force between really close particles when all smoothing is applied
 * (i.e. smoothing radius = 0). Should result in a strong repulsion.
 */
TEST_F(SmoothLennardJonesMediumDistanceTest, MediumDistance_AllSmoothing_ReallyClose) {
    SmoothLennardJonesForce force;
    force.initForce(cutoff_radius, 0);

    R3 expected_force = {-1949376.7822222223, 0.0, 0.0};
    p2.getX() = {0.5, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-0.5, 0.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}

/**
 * @brief Test the Smooth Lennard Jones Force between medium distance particles when all smoothing is applied
 * (i.e. smoothing radius = 0), when they are in each others in the potential well. Should result in no force.
 */
TEST_F(SmoothLennardJonesMediumDistanceTest, MediumDistance_AllSmoothing_PotentialWell) {
    SmoothLennardJonesForce force;
    force.initForce(cutoff_radius, 0);

    R3 expected_force = {0.0, 0.0, 0.0};
    p2.getX() = {1.1223356471450733843382863, 0.0, 0.0};
    R3 result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, expected_force, precision);

    p2.getX() = {-1.1223356471450733843382863, 0.0, 0.0};
    result = force.applyForce(p1, p2);
    EXPECT_R3_NEAR(result, -1 * expected_force, precision);
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief Tests that Newton's 3rd law holds for our force calculations.
 * (i.e. force particle 1 = -force particle 2).
 *
 * */
TEST_F(SmoothLennardJonesForceTest, Newton3) {
    p1.getEpsilon() = big_epsilon;
    p2.getEpsilon() = big_epsilon;
    p1.getSigma() = small_sigma;
    p2.getSigma() = small_sigma;
    p2.getX() = {4.0, 4.0, 0.0};
    R3 result1 = force.applyForce(p1, p2);
    R3 result2 = force.applyForce(p2, p1);
    EXPECT_R3_NEAR(result1, -1 * result2, precision);
}

/**
 * @brief Tests that the force between particles with different epsilons and sigmas gets computed correctly.
 *
 */
TEST_F(SmoothLennardJonesForceTest, MixingRules) {
    p1.getSigma() = 1.0;
    p1.getEpsilon() = 5.0;
    p2.getSigma() = 2.0;
    p2.getEpsilon() = 10.0;

    p2.getX() = {3.0, 4.0, 0.0};

    const R3 expected_force = {1.48242203050000000141839251455E-2, 1.97656040600000006135505969951E-2, 0.};
    R3 result = SmoothLennardJonesForce().applyForce(p1, p2);

    EXPECT_R3_NEAR(result, expected_force, precision);
}

}  // namespace mol_sim
