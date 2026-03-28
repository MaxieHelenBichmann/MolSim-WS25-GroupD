#include "particles/generators/CuboidGenerator.h"

#include <gtest/gtest.h>

#include <vector>

#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "testingUtils.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the CuboidGenerator.
 * Base Config/Data of the generator is as follow:
 * Position: (0,0,0)
 * Velocity: (0,0,0)
 * Num_Particles: (2,2,2)
 * mass: 1
 * distance: 1
 * average velocity: 0.1
 */
class CuboidGeneratorTest : public testing::Test {
   protected:
    SimpleContainer particle_container;
    ContainerRef particles;
    R3 position = {0.0, 0.0, 0.0};
    R3 velocity = {0.0, 0.0, 0.0};
    N3 num_particles = {2U, 2U, 2U};
    double mass = 1.0;
    double epsilon = 5.0;
    double sigma = 1.0;
    double distance = 1.0;
    double avg_velo = 0.1;
    double init_temp = 0.01;
    CuboidGenerator generator;

    CuboidGeneratorTest()
        : particle_container(),
          particles(particle_container),
          generator(position, velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma, init_temp) {}

    void SetUp() override {
        particles.clear();
        generator.generateParticles(particles, true);
    }
};
/**
 * @brief Tests that the CuboidGenerator generates the correct number of particles.
 *
 */
TEST_F(CuboidGeneratorTest, testParticleCount) { EXPECT_EQ(particles.size(), 8); }
/**
 * @brief Tests that the Particles are generated at the correct positions.
 *
 */
TEST_F(CuboidGeneratorTest, testParticlePositions) {
    std::vector<R3> expected_positions = {{0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}, {0., 1., 1.},
                                          {1., 0., 0.}, {1., 0., 1.}, {1., 1., 0.}, {1., 1., 1.}};

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
TEST_F(CuboidGeneratorTest, testParticleMass) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getM(), mass);
    }
}
/**
 * @brief Tests that the Particles have the correct sigma.
 *
 */
TEST_F(CuboidGeneratorTest, testParticleSigma) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getSigma(), sigma);
    }
}
/**
 * @brief Tests that the Particles have the correct epsilon.
 *
 */
TEST_F(CuboidGeneratorTest, testParticleEpsilon) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getEpsilon(), epsilon);
    }
}
/**
 * @brief Tests that for 2D cubes, particles have no velocity on the z-axis.
 *
 */
TEST_F(CuboidGeneratorTest, testParticleVelocity) {
    // Create a 2D cube (z=1) to test 2D Brownian motion
    SimpleContainer particle_container_2d;
    ContainerRef particles_2d(particle_container_2d);
    N3 num_particles_2d = {5U, 5U, 1U};  // Only 1 layer in z-direction makes it 2D
    CuboidGenerator generator_2d(position, velocity, num_particles_2d, {}, mass, distance, avg_velo, epsilon, sigma,
                                 init_temp);
    generator_2d.generateParticles(particles_2d, true);

    for (auto& p : particles_2d) {
        EXPECT_EQ(p.getV()[2], 0.0) << "2D cuboid should have zero z-velocity";
    }
}
/**
 * @brief Tests that the Particles velocity are correctly distributed around the initial velocity. (mean)
 *
 */
TEST_F(CuboidGeneratorTest, testVelocityDistribution) {
    // Redefine generator and particles for this test to have more particles
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N3 num_particles = {10U, 10U, 10U};
    R3 initial_velocity = {1.0, 2.0, 0.0};
    CuboidGenerator generator(position, initial_velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                              init_temp);
    generator.generateParticles(particles, true);

    R3 mean_velocity = {0.0, 0.0, 0.0};
    for (auto& p : particles) {
        mean_velocity = mean_velocity + p.getV();
    }
    mean_velocity = mean_velocity * (1.0 / static_cast<double>(particles.size()));

    // Check if mean velocity is close to initial velocity
    EXPECT_NEAR(mean_velocity[0], initial_velocity[0], 1e-1);
    EXPECT_NEAR(mean_velocity[1], initial_velocity[1], 1e-1);
    EXPECT_NEAR(mean_velocity[2], initial_velocity[2], 1e-1);
}
/**
 * @brief Tests that the Particles velocity are correctly distributed around the initial velocity and take the avg_velo
 * into account correctly. (variance)
 *
 */
TEST_F(CuboidGeneratorTest, testAverageVelocity) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N3 num_particles = {10U, 10U, 1U};  // 2D for this test
    double avg_velo = 0.5;
    double init_temp = 0.25;
    R3 initial_velocity = {1.0, 2.0, 0.0};

    CuboidGenerator generator(position, initial_velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                              init_temp);
    generator.generateParticles(particles, true);

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
    EXPECT_NEAR(variance[0], avg_velo * avg_velo, 1e-1);
    EXPECT_NEAR(variance[1], avg_velo * avg_velo, 1e-1);
    EXPECT_NEAR(variance[2], 0.0, 1e-1) << "2D cuboid should have no z-velocity variance";
}

/**
 * @brief Tests that no particles are generated when one dimension is 0
 *
 */
TEST_F(CuboidGeneratorTest, testZeroParticleGeneration) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N3 num_particles = {10U, 0U, 10U};
    CuboidGenerator generator(position, velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                              init_temp);
    generator.generateParticles(particles, true);
    EXPECT_EQ(particles.size(), 0);
}

/**
 * @brief Tests that Brownian motion can be disabled and particles have exact initial velocity
 *
 */
TEST_F(CuboidGeneratorTest, testBrownianMotionDisabled) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N3 num_particles = {5U, 5U, 5U};
    R3 initial_velocity = {1.0, 2.0, 3.0};
    CuboidGenerator generator(position, initial_velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                              init_temp);

    // Generate particles WITHOUT Brownian motion
    generator.generateParticles(particles, false);

    // All particles should have exactly the initial velocity (no randomness)
    for (auto& p : particles) {
        EXPECT_R3_EQ(p.getV(), initial_velocity);
    }
}

/**
 * @brief Tests that thermostat mode (use_init_temp) correctly calculates velocity from temperature
 *
 */
TEST_F(CuboidGeneratorTest, testThermostatMode) {
    SimpleContainer particle_container_temp;
    SimpleContainer particle_container_avg;
    ContainerRef particles_temp(particle_container_temp);
    ContainerRef particles_avg(particle_container_avg);
    N3 num_particles = {10U, 10U, 10U};  // 3D cuboid
    R3 initial_velocity = {0.0, 0.0, 0.0};

    double init_temp_val = 100.0;  // Higher temperature
    double avg_velo_val = 0.1;     // Low average velocity
    double mass_val = 1.0;

    CuboidGenerator generator(position, initial_velocity, num_particles, {}, mass_val, distance, avg_velo_val, epsilon,
                              sigma, init_temp_val);

    // Generate with use_init_temp = true (should use temperature)
    generator.generateParticles(particles_temp, true, true);

    // Generate with use_init_temp = false (should use avg_velo)
    generator.generateParticles(particles_avg, true, false);

    // Calculate average velocity magnitudes
    double avg_vel_magnitude_temp = 0.0;
    double avg_vel_magnitude_avg = 0.0;

    for (auto& p : particles_temp) {
        R3 v = p.getV();
        avg_vel_magnitude_temp += v.euclidNorm();
    }
    avg_vel_magnitude_temp /= static_cast<double>(particles_temp.size());

    for (auto& p : particles_avg) {
        R3 v = p.getV();
        avg_vel_magnitude_avg += v.euclidNorm();
    }
    avg_vel_magnitude_avg /= static_cast<double>(particles_avg.size());

    // For 3D Maxwell-Boltzmann distribution with parameter avg_v:
    // Each component ~ N(0, avg_v²)
    // Expected magnitude = avg_v * sqrt(8/pi) ≈ avg_v * 1.5958

    const double avg_v_from_temp = sqrt(init_temp_val / mass_val);              // = 10
    const double expected_magnitude_temp = avg_v_from_temp * sqrt(8.0 / M_PI);  // ≈ 15.96
    const double expected_magnitude_avg = avg_velo_val * sqrt(8.0 / M_PI);      // ≈ 0.16

    // When using temperature, velocity magnitude should be ~sqrt(T/m) * sqrt(8/pi)
    // When using avg_velo, velocity magnitude should be ~avg_velo * sqrt(8/pi)
    EXPECT_NEAR(avg_vel_magnitude_temp, expected_magnitude_temp, 0.5)
        << "Thermostat mode should calculate velocity from temperature";
    EXPECT_NEAR(avg_vel_magnitude_avg, expected_magnitude_avg, 0.02)
        << "Non-thermostat mode should use avg_velo parameter";
}

/**
 * @brief Tests that Brownian motion disabled with thermostat still uses exact initial velocity
 *
 */
TEST_F(CuboidGeneratorTest, testBrownianDisabledWithThermostat) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N3 num_particles = {5U, 5U, 5U};
    R3 initial_velocity = {1.5, 2.5, 3.5};
    double init_temp_val = 100.0;

    CuboidGenerator generator(position, initial_velocity, num_particles, {}, mass, distance, avg_velo, epsilon, sigma,
                              init_temp_val);

    // Generate with Brownian motion disabled but thermostat enabled
    generator.generateParticles(particles, false, true);

    // Even with thermostat, no Brownian motion means exact initial velocity
    for (auto& p : particles) {
        EXPECT_R3_EQ(p.getV(), initial_velocity)
            << "Disabled Brownian motion should result in exact initial velocity regardless of thermostat";
    }
}

}  // namespace mol_sim
