#include "particles/generators/DiscGenerator.h"

#include <gtest/gtest.h>

#include <vector>

#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "testingUtils.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the DiskGenerator.
 * Base Config/Data of the generator is as follow:
 * Position: (0,0,0)
 * Velocity: (0,0,0)
 * radius: 2
 * mass: 1
 * distance: 1
 * average velocity: 0.1
 */

class DiscGeneratorTest : public testing::Test {
   protected:
    SimpleContainer particle_container;
    ContainerRef particles;
    R3 position = {0.0, 0.0, 0.0};
    R3 velocity = {0.0, 0.0, 0.0};
    size_t radius = 2;
    double mass = 1.0;
    double epsilon = 5.0;
    double sigma = 1.0;
    double distance = 1.0;
    double avg_velo = 0.1;
    double init_temp = 0.01;
    double precision = 1e-1;
    DiscGenerator generator;

    DiscGeneratorTest()
        : particle_container(),
          particles(particle_container),
          generator(position, velocity, radius, mass, distance, avg_velo, epsilon, sigma, init_temp) {}

    void SetUp() override {
        particles.clear();
        generator.generateParticles(particles);
    }
    /*
     *
     ***
     *****
     ***
     *
     */
};
/**
 * @brief Tests that the DiscGenerator generates the correct number of particles.
 *
 */
TEST_F(DiscGeneratorTest, testParticleCount) { EXPECT_EQ(particles.size(), 13); }
/**
 * @brief Tests that the Particles are generated at the correct positions.
 *
 */
TEST_F(DiscGeneratorTest, testParticlePositions) {
    std::vector<R3> expected_positions = {
        {0., 0., 0.},  {-1., 0., 0.},  {1., 0., 0.},  {0., -1., 0.}, {0., 1., 0.},  {1., 1., 0.}, {-1., 1., 0.},
        {1., -1., 0.}, {-1., -1., 0.}, {-2., 0., 0.}, {2., 0., 0.},  {0., -2., 0.}, {0., 2., 0.},
    };

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
 *
 */
TEST_F(DiscGeneratorTest, testParticleMass) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getM(), mass);
    }
}
/**
 * @brief Tests that the Particles have the correct sigma.
 *
 */
TEST_F(DiscGeneratorTest, testParticleSigma) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getSigma(), sigma);
    }
}
/**
 * @brief Tests that the Particles have the correct epsilon.
 *
 */
TEST_F(DiscGeneratorTest, testParticleEpsilon) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getEpsilon(), epsilon);
    }
}
/**
 * @brief Tests that the Particles have no velocity on the z-axis.
 *
 */
TEST_F(DiscGeneratorTest, testParticleVelocity) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getV()[2], 0.0);
    }
}
/**
 * @brief Tests that the Particles velocity are correctly distributed around the initial velocity. (mean)
 *
 */
TEST_F(DiscGeneratorTest, testVelocityDistribution) {
    // Redefine generator and particles for this test to have more particles
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    radius = 10;
    R3 initial_velocity = {1.0, 2.0, 0.0};
    DiscGenerator generator(position, initial_velocity, radius, mass, distance, avg_velo, epsilon, sigma, init_temp);
    generator.generateParticles(particles);

    R3 mean_velocity = {0.0, 0.0, 0.0};
    for (auto& p : particles) {
        mean_velocity = mean_velocity + p.getV();
    }
    mean_velocity = mean_velocity * (1.0 / static_cast<double>(particles.size()));

    // Check if mean velocity is close to initial velocity
    EXPECT_R3_NEAR(mean_velocity, initial_velocity, precision);
}
/**
 * @brief Tests that the Particles velocity are correctly distributed around the initial velocity and take the avg_velo
 * into account correctly. (variance)
 *
 */
TEST_F(DiscGeneratorTest, testAverageVelocity) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    radius = 10;
    double avg_velo = 0.5;
    R3 initial_velocity = {1.0, 2.0, 0.0};

    DiscGenerator generator(position, initial_velocity, radius, mass, distance, avg_velo, epsilon, sigma, init_temp);
    generator.generateParticles(particles);

    // Calculate mean velocity
    R3 mean_velocity = {0.0, 0.0, 0.0};
    for (auto& p : particles) {
        mean_velocity = mean_velocity + p.getV();
    }
    mean_velocity = mean_velocity * (1.0 / static_cast<double>(particles.size()));
    // Calculate variance
    R3 variance = {0.0, 0.0, 0.0};
    for (auto& p : particles) {
        R3 v = p.getV();
        variance[0] += (v[0] - mean_velocity[0]) * (v[0] - mean_velocity[0]);
        variance[1] += (v[1] - mean_velocity[1]) * (v[1] - mean_velocity[1]);
        variance[2] += (v[2] - mean_velocity[2]) * (v[2] - mean_velocity[2]);
    }
    variance = variance * (1.0 / static_cast<double>(particles.size()));

    // Check if variance is close to avg_velo^2 for 2D, and 0 for the 3rd dimension
    R3 expected_var = {avg_velo * avg_velo, avg_velo * avg_velo, 0.};
    EXPECT_R3_NEAR(variance, expected_var, precision);
}

/**
 * @brief Tests that no particles are generated when radius is 0
 *
 */
TEST_F(DiscGeneratorTest, testZeroParticleGeneration) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    radius = 0;
    DiscGenerator generator(position, velocity, radius, mass, distance, avg_velo, epsilon, sigma, init_temp);
    generator.generateParticles(particles);
    EXPECT_EQ(particles.size(), 0);
}

}  // namespace mol_sim
