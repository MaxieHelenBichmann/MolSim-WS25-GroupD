#include "particles/generators/MembraneGenerator.h"

#include <gtest/gtest.h>

#include <vector>

#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the MembraneGenerator.
 * Base Config/Data of the generator is as follows:
 * Position: (0,0,0)
 * Velocity: (0,0,0)
 * Num_Particles: (3,3) - 2D membrane
 * mass: 1
 * distance: 1
 * average velocity: 0.1
 */
class MembraneGeneratorTest : public testing::Test {
   protected:
    SimpleContainer particle_container;
    ContainerRef particles;
    R3 position = {0.0, 0.0, 0.0};
    R3 velocity = {0.0, 0.0, 0.0};
    N2 num_particles = {3U, 3U};
    double mass = 1.0;
    double epsilon = 5.0;
    double sigma = 1.0;
    double distance = 1.0;
    double avg_velo = 0.1;
    double init_temp = 0.01;
    MembraneGenerator generator;

    MembraneGeneratorTest()
        : particle_container(),
          particles(particle_container),
          generator(position, velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma, init_temp) {}

    void SetUp() override {
        particles.clear();
        generator.generateParticles(particles);
    }
};

/**
 * @brief Tests that the MembraneGenerator generates the correct number of particles.
 */
TEST_F(MembraneGeneratorTest, testParticleCount) { EXPECT_EQ(particles.size(), 9); }

/**
 * @brief Tests that the Particles are generated at the correct positions.
 */
TEST_F(MembraneGeneratorTest, testParticlePositions) {
    std::vector<R3> expected_positions = {{0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}, {0., 1., 0.}, {1., 1., 0.},
                                          {2., 1., 0.}, {0., 2., 0.}, {1., 2., 0.}, {2., 2., 0.}};

    for (R3 pos : expected_positions) {
        bool found = false;
        for (auto& p : particles) {
            if (p.getX() == pos) {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "Particle at position " << pos << " not found.";
    }
}

/**
 * @brief Tests that the Particles have the correct mass.
 */
TEST_F(MembraneGeneratorTest, testParticleMass) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getM(), mass);
    }
}

/**
 * @brief Tests that the Particles have the correct sigma.
 */
TEST_F(MembraneGeneratorTest, testParticleSigma) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getSigma(), sigma);
    }
}

/**
 * @brief Tests that the Particles have the correct epsilon.
 */
TEST_F(MembraneGeneratorTest, testParticleEpsilon) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getEpsilon(), epsilon);
    }
}

/**
 * @brief Tests that membrane particles have the correct default type (2).
 */
TEST_F(MembraneGeneratorTest, testParticleType) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getType(), 2);
    }
}

/**
 * @brief Tests that target particles are marked with type 4.
 */
TEST_F(MembraneGeneratorTest, testTargetParticleType) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    std::vector<N2> targets = {{1, 1}, {2, 2}};
    MembraneGenerator generator(position, velocity, num_particles, targets, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);

    int target_count = 0;
    int normal_count = 0;
    for (auto& p : particles) {
        if (p.getType() == 4) {
            target_count++;
        } else if (p.getType() == 2) {
            normal_count++;
        }
    }
    EXPECT_EQ(target_count, 2);
    EXPECT_EQ(normal_count, 7);
}

/**
 * @brief Tests that the CORRECT particles at specified coordinates are marked as targets.
 * Validates that particles at target coordinates have type 4 and correct positions.
 */
TEST_F(MembraneGeneratorTest, testCorrectParticlesAreTargets) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N2 num_particles = {5U, 5U};
    std::vector<N2> targets = {{0, 0}, {2, 2}, {4, 4}};
    R3 position = {10.0, 20.0, 0.0};
    double distance = 1.5;

    MembraneGenerator generator(position, velocity, num_particles, targets, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);

    // Expected target positions in world coordinates
    std::vector<R3> expected_target_positions = {
        {10.0, 20.0, 0.0},  // (0,0)
        {13.0, 23.0, 0.0},  // (2,2)
        {16.0, 26.0, 0.0}   // (4,4)
    };

    int targets_found = 0;
    for (auto& p : particles) {
        if (p.getType() == 4) {
            targets_found++;
            // Verify this particle is at one of the expected target positions
            bool position_matches = false;
            for (const auto& expected_pos : expected_target_positions) {
                if (p.getX() == expected_pos) {
                    position_matches = true;
                    break;
                }
            }
            EXPECT_TRUE(position_matches) << "Target particle found at unexpected position: " << p.getX();
        }
    }
    EXPECT_EQ(targets_found, 3) << "Not all target particles were found with type 4";

    // Verify non-target particles have type 2
    for (auto& p : particles) {
        if (p.getType() != 4) {
            EXPECT_EQ(p.getType(), 2) << "Non-target particle has incorrect type at position: " << p.getX();
        }
    }
}

/**
 * @brief Tests that all target coordinates result in target particles, even at boundaries.
 */
TEST_F(MembraneGeneratorTest, testTargetParticlesAtBoundaries) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N2 num_particles = {4U, 4U};

    // Target all four corners
    std::vector<N2> targets = {{0, 0}, {3, 0}, {0, 3}, {3, 3}};

    MembraneGenerator generator(position, velocity, num_particles, targets, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);

    EXPECT_EQ(particles.size(), 16);

    // Count and validate corner particles are targets
    std::vector<size_t> corner_indices = {0, 3, 12, 15};  // Bottom-left, bottom-right, top-left, top-right
    for (size_t idx : corner_indices) {
        EXPECT_EQ(particles[idx].getType(), 4) << "Corner particle at index " << idx << " should be a target";
    }

    // Verify we have exactly 4 targets
    int target_count = 0;
    for (auto& p : particles) {
        if (p.getType() == 4) {
            target_count++;
        }
    }
    EXPECT_EQ(target_count, 4);
}

/**
 * @brief Tests that duplicate target coordinates don't cause issues.
 */
TEST_F(MembraneGeneratorTest, testDuplicateTargetCoordinates) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);

    // Same coordinate listed multiple times
    std::vector<N2> targets = {{1, 1}, {1, 1}, {1, 1}, {2, 2}};

    MembraneGenerator generator(position, velocity, num_particles, targets, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);

    // Should still only have 2 unique target particles
    int target_count = 0;
    for (auto& p : particles) {
        if (p.getType() == 4) {
            target_count++;
        }
    }
    EXPECT_EQ(target_count, 2) << "Duplicate coordinates should not create extra target particles";
}

/**
 * @brief Tests that the Particles have no velocity on the z-axis.
 */
TEST_F(MembraneGeneratorTest, testParticleVelocity) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getV()[2], 0.0);
    }
}

/**
 * @brief Tests that the Particles velocity are correctly distributed around the initial velocity.
 */
TEST_F(MembraneGeneratorTest, testVelocityDistribution) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N2 num_particles = {20U, 20U};
    R3 initial_velocity = {1.0, 2.0, 0.0};
    MembraneGenerator generator(position, initial_velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);

    R3 mean_velocity = {0.0, 0.0, 0.0};
    for (auto& p : particles) {
        mean_velocity = mean_velocity + p.getV();
    }
    mean_velocity = mean_velocity * (1.0 / static_cast<double>(particles.size()));

    EXPECT_NEAR(mean_velocity[0], initial_velocity[0], 1e-1);
    EXPECT_NEAR(mean_velocity[1], initial_velocity[1], 1e-1);
    EXPECT_NEAR(mean_velocity[2], initial_velocity[2], 1e-1);
}

/**
 * @brief Tests that no particles are generated when one dimension is 0.
 */
TEST_F(MembraneGeneratorTest, testZeroParticleGeneration) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N2 num_particles = {10U, 0U};
    MembraneGenerator generator(position, velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);
    EXPECT_EQ(particles.size(), 0);
}

/**
 * @brief Tests neighbor setup for a 3x3 membrane grid.
 * Verifies that corner, edge, and interior particles have correct neighbor counts.
 */
TEST_F(MembraneGeneratorTest, testNeighborCounts) {
    // For a 3x3 grid:
    // Corner particles should have 3 neighbors
    // Edge particles should have 5 neighbors
    // Center particle should have 8 neighbors

    auto count_neighbors = [](const Particle& p) {
        int count = 0;
        for (const auto& neighbor : p.getNeighbors()) {
            if (neighbor.has_value()) {
                count++;
            }
        }
        return count;
    };

    // Bottom-left corner (0,0)
    EXPECT_EQ(count_neighbors(particles[0]), 3);
    // Bottom-right corner (2,0)
    EXPECT_EQ(count_neighbors(particles[2]), 3);
    // Top-left corner (0,2)
    EXPECT_EQ(count_neighbors(particles[6]), 3);
    // Top-right corner (2,2)
    EXPECT_EQ(count_neighbors(particles[8]), 3);

    // Edge particles (non-corner)
    // Bottom edge (1,0)
    EXPECT_EQ(count_neighbors(particles[1]), 5);
    // Left edge (0,1)
    EXPECT_EQ(count_neighbors(particles[3]), 5);
    // Right edge (2,1)
    EXPECT_EQ(count_neighbors(particles[5]), 5);

    // Top edge (1,2)
    EXPECT_EQ(count_neighbors(particles[7]), 5);

    // Center particle (1,1)
    EXPECT_EQ(count_neighbors(particles[4]), 8);
}

/**
 * @brief Tests specific neighbor relationships for the center particle.
 * Verifies all 8 neighbors (4 cardinal + 4 diagonal) are correctly set.
 *
 * Grid layout for 3x3 (indices increase with y):
 *   6  7  8    (y=2, top row in physical space)
 *   3  4  5    (y=1, middle row)
 *   0  1  2    (y=0, bottom row in physical space)
 */
TEST_F(MembraneGeneratorTest, testCenterParticleNeighbors) {
    // Center particle is at index 4 (position 1,1)
    Particle& center = particles[4];
    auto& neighbors = center.getNeighbors();

    // Left neighbor (index 0): particle at (0,1) - index 3
    EXPECT_EQ(neighbors[0], 3);
    // Right neighbor (index 1): particle at (2,1) - index 5
    EXPECT_EQ(neighbors[1], 5);
    // Bottom neighbor (index 2): particle at (1,0) - index 1 (lower y)
    EXPECT_EQ(neighbors[2], 1);
    // Top neighbor (index 3): particle at (1,2) - index 7 (higher y)
    EXPECT_EQ(neighbors[3], 7);
    // Bottom-left diagonal (index 4): particle at (0,0) - index 0 (lower y, lower x)
    EXPECT_EQ(neighbors[4], 0);
    // Bottom-right diagonal (index 5): particle at (2,0) - index 2 (lower y, higher x)
    EXPECT_EQ(neighbors[5], 2);
    // Top-left diagonal (index 6): particle at (0,2) - index 6 (higher y, lower x)
    EXPECT_EQ(neighbors[6], 6);
    // Top-right diagonal (index 7): particle at (2,2) - index 8 (higher y, higher x)
    EXPECT_EQ(neighbors[7], 8);
}

/**
 * @brief Tests that neighbor relationships are bidirectional.
 * If A is a neighbor of B, then B should be a neighbor of A in the opposite direction.
 */
TEST_F(MembraneGeneratorTest, testBidirectionalNeighbors) {
    // Test left-right relationship between particles 0 and 1
    Particle& p0 = particles[0];  // (0,0)
    Particle& p1 = particles[1];  // (1,0)

    // p1 should be the right neighbor (index 1) of p0
    EXPECT_EQ(p0.getNeighbors()[1], 1);
    // p0 should be the left neighbor (index 0) of p1
    EXPECT_EQ(p1.getNeighbors()[0], 0);

    // Test bottom-top relationship between particles 0 (y=0) and 3 (y=1)
    Particle& p3 = particles[3];  // (0,1)

    // p3 should be the top neighbor (index 3) of p0
    EXPECT_EQ(p0.getNeighbors()[3], 3);
    // p0 should be the bottom neighbor (index 2) of p3
    EXPECT_EQ(p3.getNeighbors()[2], 0);
}

/**
 * @brief Tests that corner particles have null neighbors in unavailable directions.
 */
TEST_F(MembraneGeneratorTest, testCornerParticleNullNeighbors) {
    Particle& bottom_left = particles[0];  // (0,0)

    // Bottom-left corner should have no left, bottom, or bottom-diagonal neighbors
    EXPECT_EQ(bottom_left.getNeighbors()[0], std::nullopt);  // Left
    EXPECT_EQ(bottom_left.getNeighbors()[2], std::nullopt);  // Bottom (no row below)
    EXPECT_EQ(bottom_left.getNeighbors()[4], std::nullopt);  // Bottom-left diagonal
    EXPECT_EQ(bottom_left.getNeighbors()[5], std::nullopt);  // Bottom-right diagonal

    Particle& top_right = particles[8];  // (2,2)

    // Top-right corner should have no right, top, or top-diagonal neighbors
    EXPECT_EQ(top_right.getNeighbors()[1], std::nullopt);  // Right
    EXPECT_EQ(top_right.getNeighbors()[3], std::nullopt);  // Top (no row above)
    EXPECT_EQ(top_right.getNeighbors()[6], std::nullopt);  // Top-left diagonal
    EXPECT_EQ(top_right.getNeighbors()[7], std::nullopt);  // Top-right diagonal
}

/**
 * @brief Tests neighbor setup for a 2x2 minimal membrane.
 */
TEST_F(MembraneGeneratorTest, testMinimalMembrane) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N2 num_particles = {2U, 2U};
    MembraneGenerator generator(position, velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);

    EXPECT_EQ(particles.size(), 4);

    // All particles in 2x2 grid should have exactly 3 neighbors (corner particles)
    for (auto& p : particles) {
        int neighbor_count = 0;
        for (const auto& neighbor : p.getNeighbors()) {
            if (neighbor.has_value()) {
                neighbor_count++;
            }
        }
        EXPECT_EQ(neighbor_count, 3);
    }
}

/**
 * @brief Tests that particles added to an existing container maintain correct neighbor indices.
 */
TEST_F(MembraneGeneratorTest, testNeighborsWithExistingParticles) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);

    // Add some particles first
    particles.addParticle({-10., -10., 0.}, {0., 0., 0.}, 1.0, 1.0, 1.0, 0);
    particles.addParticle({-20., -20., 0.}, {0., 0., 0.}, 1.0, 1.0, 1.0, 0);

    size_t base = particles.size();
    EXPECT_EQ(base, 2);

    N2 num_particles = {2U, 2U};
    MembraneGenerator generator(position, velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                                init_temp);
    generator.generateParticles(particles);

    EXPECT_EQ(particles.size(), 6);

    // Verify that membrane particles (indices 2-5) have correct neighbor relationships
    // and don't accidentally reference the pre-existing particles (indices 0-1)
    for (size_t i = base; i < particles.size(); i++) {
        for (const auto& neighbor : particles[i].getNeighbors()) {
            if (neighbor.has_value()) {
                // All neighbors should be within the membrane particle range [base, size)
                EXPECT_GE(neighbor.value(), base) << "Membrane particle has neighbor outside membrane range";
                EXPECT_LT(neighbor.value(), particles.size()) << "Membrane particle has invalid neighbor index";
            }
        }
    }
}

}  // namespace mol_sim
