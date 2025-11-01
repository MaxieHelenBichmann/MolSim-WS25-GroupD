#include "particles/generators/CuboidGenerator.h"

#include <gtest/gtest.h>

#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"

namespace mol_sim {

class CuboidGeneratorTest : public testing::Test {
   protected:
    SimpleContainer particle_container;
    ContainerRef particles;
    R3 position = {0.0, 0.0, 0.0};
    R3 velocity = {0.0, 0.0, 0.0};
    N3 num_particles = {2U, 2U, 2U};
    double mass = 1.0;
    double distance = 1.0;
    double avg_velo = 0.1;
    CuboidGenerator generator;

    CuboidGeneratorTest()
        : particle_container(),
          particles(particle_container),
          generator(position, velocity, num_particles, mass, distance, avg_velo) {}

    void SetUp() override {
        particles.clear();
        generator.generateParticles(particles);
    }
};

TEST_F(CuboidGeneratorTest, testParticleCount) { EXPECT_EQ(particles.size(), 8); }

TEST_F(CuboidGeneratorTest, testParticlePositions) {
    // Check particle positions
    for (size_t i = 0; i < num_particles[2]; i++) {
        for (size_t j = 0; j < num_particles[1]; j++) {
            for (size_t k = 0; k < num_particles[0]; k++) {
                R3 expected_pos = {k * distance, j * distance, i * distance};
                bool found = false;
                for (auto& p : particles) {
                    if (p.getX() == expected_pos) {
                        found = true;
                        break;
                    }
                }
                EXPECT_TRUE(found) << "Particle at position " << expected_pos << " not found.";
            }
        }
    }
}

TEST_F(CuboidGeneratorTest, testParticleMass) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getM(), mass);
    }
}

TEST_F(CuboidGeneratorTest, testParticleVelocity) {
    for (auto& p : particles) {
        EXPECT_EQ(p.getV()[2], 0.0);
    }
}

TEST_F(CuboidGeneratorTest, testVelocityDistribution) {
    // Redefine generator and particles for this test to have more particles
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N3 num_particles = {10U, 10U, 10U};
    R3 initial_velocity = {1.0, 2.0, 0.0};
    CuboidGenerator generator(position, initial_velocity, num_particles, mass, distance, avg_velo);
    generator.generateParticles(particles);

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

TEST_F(CuboidGeneratorTest, testAverageVelocity) {
    SimpleContainer particle_container;
    ContainerRef particles(particle_container);
    N3 num_particles = {10U, 10U, 10U};
    double avg_velo = 0.5;
    R3 initial_velocity = {1.0, 2.0, 0.0};

    CuboidGenerator generator(position, initial_velocity, num_particles, mass, distance, avg_velo);
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
    EXPECT_NEAR(variance[0], avg_velo * avg_velo, 1e-1);
    EXPECT_NEAR(variance[1], avg_velo * avg_velo, 1e-1);
    EXPECT_NEAR(variance[2], 0.0, 1e-1);
}

}  // namespace mol_sim
