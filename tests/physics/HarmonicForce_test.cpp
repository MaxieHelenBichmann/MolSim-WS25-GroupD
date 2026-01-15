#include "physics/singleforces/HarmonicForce.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "particles/Particle.h"
#include "testingUtils.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the Harmonic Force
 * Tests the spring-like forces between membrane particle neighbors.
 */
class HarmonicForceTest : public testing::Test {
   protected:
    const double K = 300.0;  // Stiffness constant
    const double R_0 = 2.2;  // Equilibrium distance for direct neighbors
    double precision = 1e-9;
    HarmonicForce force;

    HarmonicForceTest() : force(K, R_0) {}
};

/**
 * @brief Tests that non-membrane particles (type 0) experience no harmonic force.
 */
TEST_F(HarmonicForceTest, NonMembraneParticleNoForce) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 0);  // type 0 = non-membrane

    R3 result = force.applyForce(p1);
    EXPECT_R3_EQ((R3{0.0, 0.0, 0.0}), result);
}

/**
 * @brief Tests that membrane particles with no neighbors experience no force.
 */
TEST_F(HarmonicForceTest, NoNeighborsNoForce) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // type 2 = membrane

    // All neighbors are nullptr by default
    R3 result = force.applyForce(p1);
    EXPECT_R3_EQ((R3{0.0, 0.0, 0.0}), result);
}

/**
 * @brief Tests harmonic force when particle is at equilibrium distance from one neighbor.
 */
TEST_F(HarmonicForceTest, SingleNeighborAtEquilibrium) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p2(R3{R_0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    p1.getNeighbors()[1] = &p2;

    R3 result = force.applyForce(p1);
    EXPECT_R3_NEAR((R3{0.0, 0.0, 0.0}), result, precision);
}

/**
 * @brief Tests harmonic force when particle is stretched from neighbor (distance > R_0).
 */
TEST_F(HarmonicForceTest, SingleNeighborStretched) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p2(R3{R_0 + 1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    p1.getNeighbors()[1] = &p2;

    double dist = R_0 + 1.0;
    double expected_magnitude = K * 0.5 * (dist - R_0);
    R3 expected_force = {expected_magnitude, 0.0, 0.0};

    R3 result = force.applyForce(p1);
    EXPECT_R3_NEAR(result, expected_force, precision);
}

/**
 * @brief Tests harmonic force when particle is compressed toward neighbor (distance < R_0).
 */
TEST_F(HarmonicForceTest, SingleNeighborCompressed) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p2(R3{R_0 - 0.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    p1.getNeighbors()[1] = &p2;

    double dist = R_0 - 0.5;
    double expected_magnitude = K * 0.5 * (dist - R_0);
    R3 expected_force = {expected_magnitude, 0.0, 0.0};

    R3 result = force.applyForce(p1);
    EXPECT_R3_NEAR(result, expected_force, precision);
}

/**
 * @brief Tests harmonic force with diagonal neighbor at equilibrium.
 */
TEST_F(HarmonicForceTest, DiagonalNeighborAtEquilibrium) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Diagonal equilibrium distance is sqrt(2) * R_0
    double diag_dist = std::numbers::sqrt2 * R_0;
    Particle p2(R3{diag_dist / std::numbers::sqrt2, diag_dist / std::numbers::sqrt2, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0,
                1.0, 2);

    p1.getNeighbors()[7] = &p2;  // Top-right diagonal

    R3 result = force.applyForce(p1);
    EXPECT_R3_NEAR(result, (R3{0.0, 0.0, 0.0}), precision);
}

/**
 * @brief Tests harmonic force with diagonal neighbor stretched.
 */
TEST_F(HarmonicForceTest, DiagonalNeighborStretched) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Place diagonal neighbor beyond equilibrium
    double extra = 0.5;
    double total_dist = (std::numbers::sqrt2 * R_0) + extra;
    Particle p2(R3{total_dist / std::numbers::sqrt2, total_dist / std::numbers::sqrt2, 0.0}, R3{0.0, 0.0, 0.0}, 1.0,
                1.0, 1.0, 2);

    p1.getNeighbors()[7] = &p2;  // Top-right diagonal

    double dist = total_dist;
    double expected_magnitude = K * 0.5 * (dist - std::numbers::sqrt2 * R_0) / dist;
    R3 direction = p2.getX() - p1.getX();
    R3 expected_force = expected_magnitude * direction;

    R3 result = force.applyForce(p1);
    EXPECT_R3_NEAR(result, expected_force, precision);
}

/**
 * @brief Tests that forces from multiple neighbors sum correctly.
 */
TEST_F(HarmonicForceTest, MultipleNeighborsSumForces) {
    Particle p_center(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Create four direct neighbors, all stretched by same amount
    double stretch = 0.3;
    Particle p_left(R3{-(R_0 + stretch), 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_right(R3{R_0 + stretch, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_bottom(R3{0.0, -(R_0 + stretch), 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_top(R3{0.0, R_0 + stretch, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    p_center.getNeighbors()[0] = &p_left;
    p_center.getNeighbors()[1] = &p_right;
    p_center.getNeighbors()[2] = &p_bottom;
    p_center.getNeighbors()[3] = &p_top;

    R3 result = force.applyForce(p_center);

    // Each neighbor pulls with equal force, but in opposite directions
    // Net force should be zero (symmetric configuration)
    EXPECT_R3_NEAR(result, (R3{0.0, 0.0, 0.0}), precision);
}

/**
 * @brief Tests asymmetric force when only one side is stretched.
 */
TEST_F(HarmonicForceTest, AsymmetricNeighbors) {
    Particle p_center(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Left neighbor at equilibrium, right neighbor stretched
    Particle p_left(R3{-R_0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_right(R3{R_0 + 1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    p_center.getNeighbors()[0] = &p_left;
    p_center.getNeighbors()[1] = &p_right;

    R3 result = force.applyForce(p_center);
    EXPECT_GT(result[0], 0.0);
    EXPECT_NEAR(result[1], 0.0, precision);
    EXPECT_NEAR(result[2], 0.0, precision);
}

/**
 * @brief Tests that target particles (type 4) also experience harmonic forces.
 */
TEST_F(HarmonicForceTest, TargetParticleExperiencesForce) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 4);  // type 4 = target membrane
    Particle p2(R3{R_0 + 0.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    p1.getNeighbors()[1] = &p2;

    R3 result = force.applyForce(p1);

    // Should experience attractive force
    EXPECT_GT(result[0], 0.0);
}

/**
 * @brief Tests force calculation with all 8 neighbors present.
 */
TEST_F(HarmonicForceTest, AllEightNeighbors) {
    Particle p_center(R3{5.0, 5.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Create 8 neighbors at equilibrium distances
    Particle p_left(R3{5.0 - R_0, 5.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_right(R3{5.0 + R_0, 5.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_bottom(R3{5.0, 5.0 - R_0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_top(R3{5.0, 5.0 + R_0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    double d_eq = std::numbers::sqrt2 * R_0;
    Particle p_bl(R3{5.0 - (d_eq / std::numbers::sqrt2), 5.0 - (d_eq / std::numbers::sqrt2), 0.0}, R3{0.0, 0.0, 0.0},
                  1.0, 1.0, 1.0, 2);
    Particle p_br(R3{5.0 + (d_eq / std::numbers::sqrt2), 5.0 - (d_eq / std::numbers::sqrt2), 0.0}, R3{0.0, 0.0, 0.0},
                  1.0, 1.0, 1.0, 2);
    Particle p_tl(R3{5.0 - (d_eq / std::numbers::sqrt2), 5.0 + (d_eq / std::numbers::sqrt2), 0.0}, R3{0.0, 0.0, 0.0},
                  1.0, 1.0, 1.0, 2);
    Particle p_tr(R3{5.0 + (d_eq / std::numbers::sqrt2), 5.0 + (d_eq / std::numbers::sqrt2), 0.0}, R3{0.0, 0.0, 0.0},
                  1.0, 1.0, 1.0, 2);

    p_center.getNeighbors()[0] = &p_left;
    p_center.getNeighbors()[1] = &p_right;
    p_center.getNeighbors()[2] = &p_bottom;
    p_center.getNeighbors()[3] = &p_top;
    p_center.getNeighbors()[4] = &p_bl;
    p_center.getNeighbors()[5] = &p_br;
    p_center.getNeighbors()[6] = &p_tl;
    p_center.getNeighbors()[7] = &p_tr;

    R3 result = force.applyForce(p_center);

    // All at equilibrium, symmetric configuration -> net force should be zero
    EXPECT_R3_NEAR((R3{0.0, 0.0, 0.0}), result, precision);
}

/**
 * @brief Tests that only the first 4 indices are used for direct neighbors.
 */
TEST_F(HarmonicForceTest, DirectNeighborIndicesCorrect) {
    Particle p_center(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p_right(R3{R_0 + 1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Test each direct neighbor index (0-3)
    for (size_t i = 0; i < 4; i++) {
        Particle p_test(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
        p_test.getNeighbors()[i] = &p_right;

        R3 result = force.applyForce(p_test);
        // Should have non-zero force (using direct neighbor equilibrium R_0)
        EXPECT_GT(result.euclidNorm(), 0.0) << "Direct neighbor at index " << i << " should contribute force";
    }
}

/**
 * @brief Tests that indices 4-7 are used for diagonal neighbors with different equilibrium.
 */
TEST_F(HarmonicForceTest, DiagonalNeighborIndicesCorrect) {
    Particle p_center(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Place a particle at direct neighbor equilibrium distance
    Particle p_neighbor(R3{R_0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    // Test each diagonal neighbor index (4-7)
    for (size_t i = 4; i < 8; i++) {
        Particle p_test(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
        p_test.getNeighbors()[i] = &p_neighbor;

        R3 result = force.applyForce(p_test);
        // Should have non-zero force (using diagonal equilibrium sqrt(2)*R_0, which is != R_0)
        EXPECT_GT(result.euclidNorm(), 0.0) << "Diagonal neighbor at index " << i << " should contribute force";
    }
}

/**
 * @brief Tests force magnitude scaling with stiffness constant K.
 */
TEST_F(HarmonicForceTest, StiffnessScaling) {
    Particle p1(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    Particle p2(R3{R_0 + 0.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    p1.getNeighbors()[1] = &p2;

    // Test with different K values
    HarmonicForce force_low(100.0, R_0);
    HarmonicForce force_high(500.0, R_0);

    R3 result_low = force_low.applyForce(p1);
    R3 result_high = force_high.applyForce(p1);

    // Higher K should produce proportionally higher force
    EXPECT_NEAR(result_high[0] / result_low[0], 5.0, 1e-6);
}

}  // namespace mol_sim
