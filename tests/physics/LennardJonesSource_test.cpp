#include <gtest/gtest.h>
#include <numbers>

#include "particles/Particle.h"
#include "physics/LennardJonesForce.h"

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

    LennardJonesForceTest() : p1{p1_x, p1_v, 1.0}, p2{p2_x, p2_v, 1.0} {}
};

/**
 * @brief Tests the computed force against hand computed force for a particle that has
 * a decent distance from the other particle. This means 
 * the particle should experience slight attraction.
 *
 * */
TEST_F(LennardJonesForceTest, FarAway_BigEpsSmallSig) {
    //hand computed values converted to most accurate double precision float
    //(https://binaryconvert.com/result_double.html)
    const R3 cardinal_x = {1.7881391299999999373561898261E-6, 0, 0};
    const R3 cardinal_y = {0, 1.7881391299999999373561898261E-6, 0};
    const R3 diagonal = {2.28880089070000009802280382942E-4, 2.28880089070000009802280382942E-4, 0};

    //Cardinals
    p2.getX() = {4.0, 0.0, 0.0};
    R3 result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_x);

    p2.getX() = {-4.0, 0.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_x);

    p2.getX() = {0.0, 4.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_y);

    p2.getX() = {0.0, -4.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_y);

    //Diagonal
    p2.getX() = {4.0, 4.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);

    p2.getX() = {-4.0, 4.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);

    p2.getX() = {-4.0, -4.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);

    p2.getX() = {4.0, -4.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);
}

/**
 * @brief Tests the computed force against hand computed force when particle is at bottom of potential
 * well. Force should be 0 there.
 *
 * */
TEST_F(LennardJonesForceTest, AtBottomOfWell) {
    //hand computed values converted to most accurate double precision float
    //(https://binaryconvert.com/result_double.html)

    /** POTENTIAL WELL AT RADIUS 2 */
    double sigma = 1.78179743599999995673499597615E0; //2^{5/6}
    double bottom_diagonal = std::numbers::sqrt2;
    R3 zero = {.0, .0, .0};

    //Cardinals
    p2.getX() = {2.0, 0.0, 0.0};
    R3 result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {-2.0, 0.0, 0.0};
    result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {0.0, 2.0, 0.0};
    result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {0.0, -2.0, 0.0};
    result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    //Diagonal
    p2.getX() = {bottom_diagonal, bottom_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {-bottom_diagonal, bottom_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {-bottom_diagonal, -bottom_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {bottom_diagonal, -bottom_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    /** POTENTIAL WELL AT RADIUS 1 */
    sigma = 8.90898718099999986641535087983E-1; //2^{-1/6}
    bottom_diagonal = 1.0/std::numbers::sqrt2;

    //Cardinals
    p2.getX() = {1.0, 0.0, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {-1.0, 0.0, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {0.0, 1.0, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {0.0, -1.0, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    //Diagonal
    p2.getX() = {bottom_diagonal, bottom_diagonal, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {-bottom_diagonal, bottom_diagonal, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {-bottom_diagonal, -bottom_diagonal, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);

    p2.getX() = {bottom_diagonal, -bottom_diagonal, 0.0};
    result = LennardJonesForce(small_epsilon, sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, zero);
}

/**
 * @brief Tests compute force against hand computed force for particle thats exactly at 
 * sigma of the other particle's potential well (i.e. at the 
 * zero crossing of the potential). We should experience a strong
 * repulsion at this point.
 *
 * */
TEST_F(LennardJonesForceTest, AtSigma) {
    R3 cardinal_x = {-120.0, 0.0, 0.0};
    R3 cardinal_y = {0.0, -120.0, 0.0};
    R3 diagonal = {-60.0 * std::numbers::sqrt2, -60.0 * std::numbers::sqrt2, 0.0};
    double sigma_diagonal = 1.0/std::numbers::sqrt2;

    //Cardinal
    p2.getX() = {small_sigma, 0.0, 0.0};
    R3 result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_x);

    p2.getX() = {-small_sigma, 0.0, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_x);

    p2.getX() = {0.0, small_sigma, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_y);

    p2.getX() = {0.0, -small_sigma, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, cardinal_y);

    //Diagonal 
    p2.getX() = {sigma_diagonal, sigma_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);

    p2.getX() = {-sigma_diagonal, sigma_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);

    p2.getX() = {-sigma_diagonal, -sigma_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);

    p2.getX() = {sigma_diagonal, -sigma_diagonal, 0.0};
    result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    EXPECT_EQ(result, diagonal);

}

/**
 * @brief Tests that Newton's 3rd law holds for our force calculations.
 * (i.e. force particle 1 = -force particle 2).
 *
 * */
TEST_F(LennardJonesForceTest, Newton3) {
    //Cardinals
    p2.getX() = {4.0, 0.0, 0.0};
    R3 result1 = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    R3 result2 = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);

    p2.getX() = {-4.0, 0.0, 0.0};
    result1 = LennardJonesForce(small_epsilon, small_sigma).calculateForce(p1, p2);
    result2 = LennardJonesForce(small_epsilon, small_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);

    p2.getX() = {0.0, 4.0, 0.0};
    result1 = LennardJonesForce(big_epsilon, big_sigma).calculateForce(p1, p2);
    result2 = LennardJonesForce(big_epsilon, big_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);

    p2.getX() = {0.0, -4.0, 0.0};
    result1 = LennardJonesForce(big_epsilon, big_sigma).calculateForce(p1, p2);
    result2 = LennardJonesForce(big_epsilon, big_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);

    //Diagonal
    p2.getX() = {4.0, 4.0, 0.0};
    result1 = LennardJonesForce(small_epsilon, big_sigma).calculateForce(p1, p2);
    result2 = LennardJonesForce(small_epsilon, big_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);

    p2.getX() = {-4.0, 4.0, 0.0};
    result1 = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    result2 = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);

    p2.getX() = {-4.0, -4.0, 0.0};
    result1 = LennardJonesForce(big_epsilon, big_sigma).calculateForce(p1, p2);
    result2 = LennardJonesForce(big_epsilon, big_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);

    p2.getX() = {4.0, -4.0, 0.0};
    result1 = LennardJonesForce(small_epsilon, small_sigma).calculateForce(p1, p2);
    result2 = LennardJonesForce(small_epsilon, small_sigma).calculateForce(p2, p1);
    EXPECT_EQ(result1, -1 * result2);
}

/**
 * @brief Tests that the force is 0 when the particles are at the same position.
 *
 */
TEST_F(LennardJonesForceTest, TestForceCalcZeroDistance) {
    p2 = {p1_x, p2_v, p2_mass};
    R3 result = LennardJonesForce(big_epsilon, small_sigma).calculateForce(p1, p2);
    R3 expected = {0.0, 0.0, 0.0};
    EXPECT_EQ(result, expected);
}

}  // namespace mol_sim
