#include "physics/singleforces/HarmonicForce.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "particles/Particle.h"
#include "particles/container/SimpleContainer.h"
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
    SimpleContainer container;

    HarmonicForceTest() : force(K, R_0), container({10.0, 10.0, 10.0}, std::numeric_limits<double>::infinity()) {
        force.setContainer(container);
    }
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
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    container.addParticle(R3{R_0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    container[0].getNeighbors()[1] = 1;

    R3 result = force.applyForce(container[0]);
    EXPECT_R3_NEAR((R3{0.0, 0.0, 0.0}), result, precision);
}

/**
 * @brief Tests harmonic force when particle is stretched from neighbor (distance > R_0).
 */
TEST_F(HarmonicForceTest, SingleNeighborStretched) {
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    container.addParticle(R3{R_0 + 1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    container[0].getNeighbors()[1] = 1;

    const R3 expected_force = {150., 0.0, 0.0};

    R3 result = force.applyForce(container[0]);
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Tests harmonic force when particle is compressed toward neighbor (distance < R_0).
 */
TEST_F(HarmonicForceTest, SingleNeighborCompressed) {
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    container.addParticle(R3{R_0 - 0.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    container[0].getNeighbors()[1] = 1;

    const R3 expected_force = {-75., 0.0, 0.0};

    R3 result = force.applyForce(container[0]);
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Tests harmonic force with diagonal neighbor at equilibrium.
 */
TEST_F(HarmonicForceTest, DiagonalNeighborAtEquilibrium) {
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    container.addParticle(R3{R_0, R_0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    container[0].getNeighbors()[7] = 1;  // Top-right diagonal

    R3 result = force.applyForce(container[0]);
    EXPECT_R3_NEAR((R3{0.0, 0.0, 0.0}), result, precision);
}

/**
 * @brief Tests harmonic force with diagonal neighbor stretched.
 */
TEST_F(HarmonicForceTest, DiagonalNeighborStretched) {
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    container.addParticle(R3{R_0 + (0.5 * std::numbers::sqrt2), R_0 + (0.5 * std::numbers::sqrt2), 0.0},
                          R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    container[0].getNeighbors()[7] = 1;  // Top-right diagonal

    // Diagonal equilibrium: sqrt(2)*R_0
    // Actual distance: sqrt(2) * (R_0 + 0.5*sqrt(2))
    // Stretch: sqrt(2) * 0.5 * sqrt(2) = 1.0
    // Force magnitude: K * 1.0 / 2 = 150
    // Direction: 45° (equal x,y components): 150/sqrt(2) ≈ 106.066 each
    const R3 expected_force = {150.0 / std::numbers::sqrt2, 150.0 / std::numbers::sqrt2, 0.0};
    R3 result = force.applyForce(container[0]);
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Tests that forces from multiple neighbors sum correctly.
 */
TEST_F(HarmonicForceTest, MultipleNeighborsSumForces) {
    double stretch = 0.3;
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);               // 0: center
    container.addParticle(R3{-(R_0 + stretch), 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 1: left
    container.addParticle(R3{R_0 + stretch, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);     // 2: right
    container.addParticle(R3{0.0, -(R_0 + stretch), 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 3: bottom
    container.addParticle(R3{0.0, R_0 + stretch, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);     // 4: top

    container[0].getNeighbors()[0] = 1;
    container[0].getNeighbors()[1] = 2;
    container[0].getNeighbors()[2] = 3;
    container[0].getNeighbors()[3] = 4;

    R3 result = force.applyForce(container[0]);

    // Each neighbor pulls with equal force, but in opposite directions
    // Net force should be zero (symmetric configuration)
    EXPECT_R3_NEAR((R3{0.0, 0.0, 0.0}), result, precision);
}

/**
 * @brief Tests asymmetric force when only one side is stretched.
 */
TEST_F(HarmonicForceTest, AsymmetricNeighbors) {
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);        // 0: center
    container.addParticle(R3{-R_0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);       // 1: left
    container.addParticle(R3{R_0 + 1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 2: right

    container[0].getNeighbors()[0] = 1;
    container[0].getNeighbors()[1] = 2;

    R3 result = force.applyForce(container[0]);
    // Left at equilibrium: 0 force, Right stretched by 1.0: K*(1.0)/2 = 150 force to the right
    const R3 expected_force = {150.0, 0.0, 0.0};
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Tests that target particles (type 4) also experience harmonic forces.
 */
TEST_F(HarmonicForceTest, TargetParticleExperiencesForce) {
    container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 4);  // type 4 = target membrane
    container.addParticle(R3{R_0 + 0.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);

    container[0].getNeighbors()[1] = 1;

    R3 result = force.applyForce(container[0]);

    // Stretched by 0.5: K * (0.5) / 2 = 300 * 0.5 / 2 = 75
    const R3 expected_force = {75.0, 0.0, 0.0};
    EXPECT_R3_NEAR(expected_force, result, precision);
}

/**
 * @brief Tests force calculation with all 8 neighbors present.
 */
TEST_F(HarmonicForceTest, AllEightNeighbors) {
    double d_eq = std::numbers::sqrt2 * R_0;
    container.addParticle(R3{5.0, 5.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);        // 0: center
    container.addParticle(R3{5.0 - R_0, 5.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 1: left
    container.addParticle(R3{5.0 + R_0, 5.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 2: right
    container.addParticle(R3{5.0, 5.0 - R_0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 3: bottom
    container.addParticle(R3{5.0, 5.0 + R_0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 4: top
    container.addParticle(R3{5.0 - (d_eq / std::numbers::sqrt2), 5.0 - (d_eq / std::numbers::sqrt2), 0.0},
                          R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 5: bl
    container.addParticle(R3{5.0 + (d_eq / std::numbers::sqrt2), 5.0 - (d_eq / std::numbers::sqrt2), 0.0},
                          R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 6: br
    container.addParticle(R3{5.0 - (d_eq / std::numbers::sqrt2), 5.0 + (d_eq / std::numbers::sqrt2), 0.0},
                          R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 7: tl
    container.addParticle(R3{5.0 + (d_eq / std::numbers::sqrt2), 5.0 + (d_eq / std::numbers::sqrt2), 0.0},
                          R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);  // 8: tr

    container[0].getNeighbors()[0] = 1;
    container[0].getNeighbors()[1] = 2;
    container[0].getNeighbors()[2] = 3;
    container[0].getNeighbors()[3] = 4;
    container[0].getNeighbors()[4] = 5;
    container[0].getNeighbors()[5] = 6;
    container[0].getNeighbors()[6] = 7;
    container[0].getNeighbors()[7] = 8;

    R3 result = force.applyForce(container[0]);

    // All at equilibrium, symmetric configuration -> net force should be zero
    EXPECT_R3_NEAR((R3{0.0, 0.0, 0.0}), result, precision);
}

/**
 * @brief Tests that only the first 4 indices are used for direct neighbors.
 */
TEST_F(HarmonicForceTest, DirectNeighborIndicesCorrect) {
    SimpleContainer test_container({10.0, 10.0, 10.0}, std::numeric_limits<double>::infinity());
    test_container.addParticle(R3{R_0 + 1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    HarmonicForce test_force(K, R_0);
    test_force.setContainer(test_container);

    // Test each direct neighbor index (0-3)
    for (size_t i = 0; i < 4; i++) {
        SimpleContainer tmp_container({10.0, 10.0, 10.0}, std::numeric_limits<double>::infinity());
        tmp_container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
        tmp_container.addParticle(R3{R_0 + 1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
        tmp_container[0].getNeighbors()[i] = 1;

        HarmonicForce tmp_force(K, R_0);
        tmp_force.setContainer(tmp_container);

        R3 result = tmp_force.applyForce(tmp_container[0]);
        // Direct neighbor stretched by 1.0: expected force magnitude = K * 1.0 / 2 = 150
        const double expected_magnitude = 150.0;
        EXPECT_NEAR(expected_magnitude, result.euclidNorm(), precision)
            << "Direct neighbor at index " << i << " should contribute correct force magnitude";
    }
}

/**
 * @brief Tests that indices 4-7 are used for diagonal neighbors with different equilibrium.
 */
TEST_F(HarmonicForceTest, DiagonalNeighborIndicesCorrect) {
    // Test each diagonal neighbor index (4-7)
    for (size_t i = 4; i < 8; i++) {
        SimpleContainer tmp_container({10.0, 10.0, 10.0}, std::numeric_limits<double>::infinity());
        tmp_container.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
        tmp_container.addParticle(R3{R_0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
        tmp_container[0].getNeighbors()[i] = 1;

        HarmonicForce tmp_force(K, R_0);
        tmp_force.setContainer(tmp_container);

        R3 result = tmp_force.applyForce(tmp_container[0]);
        // Diagonal equilibrium is sqrt(2)*R_0, but distance is R_0
        // Compression: sqrt(2)*R_0 - R_0 = R_0*(sqrt(2) - 1) ≈ 0.9113
        // Force magnitude: K * compression / 2 = 300 * 0.9113 / 2 ≈ 136.7
        const double expected_magnitude = K * R_0 * (std::numbers::sqrt2 - 1.0) / 2.0;
        EXPECT_NEAR(expected_magnitude, result.euclidNorm(), precision)
            << "Diagonal neighbor at index " << i << " should contribute correct force magnitude";
    }
}

/**
 * @brief Tests force magnitude scaling with stiffness constant K.
 */
TEST_F(HarmonicForceTest, StiffnessScaling) {
    SimpleContainer test_container_low({10.0, 10.0, 10.0}, std::numeric_limits<double>::infinity());
    test_container_low.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    test_container_low.addParticle(R3{R_0 + 0.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    test_container_low[0].getNeighbors()[1] = 1;

    SimpleContainer test_container_high({10.0, 10.0, 10.0}, std::numeric_limits<double>::infinity());
    test_container_high.addParticle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    test_container_high.addParticle(R3{R_0 + 0.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 2);
    test_container_high[0].getNeighbors()[1] = 1;

    // Test with different K values
    HarmonicForce force_low(100.0, R_0);
    HarmonicForce force_high(500.0, R_0);
    force_low.setContainer(test_container_low);
    force_high.setContainer(test_container_high);

    R3 result_low = force_low.applyForce(test_container_low[0]);
    R3 result_high = force_high.applyForce(test_container_high[0]);

    // Higher K should produce proportionally higher force
    EXPECT_NEAR(5.0, result_high[0] / result_low[0], 1e-6);
}

}  // namespace mol_sim
