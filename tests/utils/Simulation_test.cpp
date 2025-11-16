#include "utils/Simulation.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <physics/GravitationalForce.h>
#include <physics/LennardJonesForce.h>

#include "particles/Particle.h"
#include "particles/container/SimpleContainer.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Mocks ForceSource.applyForce so we use better test values.
 *
 */
class ForceMock {
   public:
    // Mocked method for applying force between two particles.
    MOCK_METHOD(R3, applyForce, (const Particle& p1, const Particle& p2), (const));
};

/**
 * @brief Test fixture for the following tests. Contains two attributes
 that will be used throughout the tests.
 1) 'particles' which is a SimpleContainer containing the particles we perform the tests on
 2) 'settings' which is a SettingsParam for setting simulation parameters
 3)
 */
class SimulationTest : public testing::Test {
   protected:
    SimpleContainer particles;
    SettingsParam settings;
    ForceMock mock;
    GravitationalForce grav_force;
    LennardJonesForce lj_force;

   public:
    void SetUp() override {
        settings.setDefaults();
        particles.clear();
    }
};

/**
 * @brief Tests basic functionality of calculateX with respect to a
 delta_t != 1 and an initial velocity.
 *
 */
TEST_F(SimulationTest, calculateX_dt) {
    settings.delta_t = 0.5;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 0.0};
    particles.addParticle(x, v, 1.0, 5., 1.);
    Simulation<SimpleContainer, GravitationalForce> simulation(particles, grav_force, settings);
    simulation.calculateX<SimpleContainer>();
    R3 expected = {5.0, -10.0, 0.0};
    EXPECT_EQ(particles[0].getX(), expected);
}

/**
 * @brief Tests calculateX when a particle has an initial velocity
 and an initial force that's acting on it.
 *
 */
TEST_F(SimulationTest, calculateX_force) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, 0.0, 0.0};
    R3 f = {0.0, 10.0, 0.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    particles.addParticle(p);
    Simulation<SimpleContainer, GravitationalForce> simulation(particles, grav_force, settings);
    simulation.calculateX<SimpleContainer>();
    R3 expected = {10.0, 5.0, 0.0};
    EXPECT_EQ(particles[0].getX(), expected);
}

/**
 * @brief Tests basic functionality of calculateV for a single particle
 when no force is acting on it.
 *
 */
TEST_F(SimulationTest, calculateV_no_force) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    particles.addParticle(x, v, 1.0, 5., 1.);
    Simulation<SimpleContainer, GravitationalForce> simulation(particles, grav_force, settings);
    simulation.calculateV();
    EXPECT_EQ(particles[0].getV(), v);
}

/**
 * @brief Tests basic functionalty of calculateV for a single particle
 when there is a force acting on it.
 *
 */
TEST_F(SimulationTest, calculateV_simple) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 f = {1.0, 1.0, 1.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    particles.addParticle(p);
    Simulation<SimpleContainer, GravitationalForce> simulation(particles, grav_force, settings);
    simulation.calculateV();
    R3 expected = {10.5, -19.5, 30.5};
    EXPECT_EQ(particles[0].getV(), expected);
}

/**
 * @brief Tests a more complex setting of calculateV where the
 particle has both a non-zero 'old_f' and 'f'.
 *
 */
TEST_F(SimulationTest, calculateV_complex) {
    settings.delta_t = 0.5;

    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 old_f = {-12.0, 15.0, -20.5};
    R3 f = {20.5, -10.0, 100.0};
    Particle p(x, v, 2.0, 5., 1.);
    p.getF() = f;
    p.getOldF() = old_f;
    particles.addParticle(p);
    Simulation<SimpleContainer, GravitationalForce> simulation(particles, grav_force, settings);
    simulation.calculateV();
    R3 expected = {11.0625, -19.375, 39.9375};
    EXPECT_EQ(particles[0].getV(), expected);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctly
 especially in regards to Newton's third law. The inter-particle force
 is 1 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_simple2_pairwise) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({.0, .0, .0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({.0, .0, .0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    particles.addParticle(p1);
    particles.addParticle(p2);
    Simulation<SimpleContainer, ForceMock> simulation(particles, mock, settings);
    EXPECT_CALL(mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    simulation.calculateF();
    EXPECT_EQ(particles[0].getF(), f12);
    EXPECT_EQ(particles[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctly
 especially in regards to Newton's third law. The inter-particle force
 is 3 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_complex2_pairwise) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({.0, .0, .0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({.0, .0, .0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {102.52, -51.3, 135.711};
    EXPECT_CALL(mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    particles.addParticle(p1);
    particles.addParticle(p2);
    Simulation<SimpleContainer, ForceMock> simulation(particles, mock, settings);
    simulation.calculateF();
    EXPECT_EQ(particles[0].getF(), f12);
    EXPECT_EQ(particles[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly.
 Each inter-particle force is 1 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_simple3_pairwise) {
    // initial positions and velocities irrelevant since we're mocking
    Particle p1({1.0, .0, .0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p2({2.0, .0, .0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p3({3.0, .0, .0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p12({1.0, .0, .0}, {.0, .0, .0}, {10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p22({2.0, .0, .0}, {.0, .0, .0}, {-10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p32({3.0, .0, .0}, {.0, .0, .0}, {-20.0, 0.0, 0.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    R3 f13 = {20.0, 0.0, 0.0};
    R3 f23 = {-10.0, 0.0, 0.0};
    EXPECT_CALL(mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    particles.addParticle(p1);
    particles.addParticle(p2);
    particles.addParticle(p3);
    Simulation<SimpleContainer, ForceMock> simulation(particles, mock, settings);
    simulation.calculateF();
    R3 expected1 = {30.0, 0.0, 0.0};
    R3 expected2 = {-20.0, 0.0, 0.0};
    R3 expected3 = {-10.0, 0.0, 0.0};
    EXPECT_EQ(particles[0].getF(), expected1);
    EXPECT_EQ(particles[1].getF(), expected2);
    EXPECT_EQ(particles[2].getF(), expected3);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly.
 Each of the inter-particle forces are 3 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_complex3_pairwise) {
    auto force = std::make_unique<ForceMock>();
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, .0, .0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p2({2.0, .0, .0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p3({3.0, .0, .0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p12({1.0, .0, .0}, {.0, .0, .0}, {10.0, 5.0, 6.0}, 1.0, 5., 1.);
    Particle p22({2.0, .0, .0}, {.0, .0, .0}, {-10.0, -5.0, -6.0}, 1.0, 5., 1.);
    Particle p32({3.0, .0, .0}, {.0, .0, .0}, {-20.0, -10.0, -7.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 5.0, 6.0};
    R3 f13 = {20.0, 10.0, 7.0};
    R3 f23 = {-10.0, -5.0, -4.0};
    EXPECT_CALL(mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    particles.addParticle(p1);
    particles.addParticle(p2);
    particles.addParticle(p3);
    Simulation<SimpleContainer, ForceMock> simulation(particles, mock, settings);

    simulation.calculateF();
    R3 expected1 = {30.0, 15.0, 13.0};
    R3 expected2 = {-20.0, -10.0, -10.0};
    R3 expected3 = {-10.0, -5.0, -3.0};
    EXPECT_EQ(particles[0].getF(), expected1);
    EXPECT_EQ(particles[1].getF(), expected2);
    EXPECT_EQ(particles[2].getF(), expected3);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 gravitational forces and 2 particles.
 *
 */
TEST_F(SimulationTest, run_gravitational_timestep) {
    settings.delta_t = 0.5;
    settings.end_time = 0.5;
    Particle p1({.0, .0, .0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({-1.0, .0, .0}, {10.0, .0, .0}, {.0, .0, .0}, 0.5, 5., 1.);
    particles.addParticle(p1);
    particles.addParticle(p2);
    Simulation<SimpleContainer, GravitationalForce> simulation(particles, grav_force, settings);
    simulation.run();
    Particle p1_expect({.0, .0, .0}, {0.0078125, .0, .0}, {0.03125, .0, .0}, 1.0, 5., 1.);
    Particle p2_expect({4.0, .0, .0}, {9.984375, .0, .0}, {-0.03125, .0, .0}, 0.5, 5., 1.);
    EXPECT_EQ(particles[0], p1_expect);
    EXPECT_EQ(particles[1], p2_expect);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 Lennard-Jones forces and 2 particles.
 *
 */
TEST_F(SimulationTest, run_lennardjones_timestep) {
    settings.delta_t = 0.5;
    settings.end_time = 0.5;
    Particle p1({.0, .0, .0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({-4.0, .0, .0}, {10.0, .0, .0}, {.0, .0, .0}, 0.5, 5., 1.);
    particles.addParticle(p1);
    particles.addParticle(p2);
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    Particle p1_expect({.0, .0, .0}, {-30, .0, .0}, {-120, .0, .0}, 1.0, 5., 1.);
    Particle p2_expect({1.0, .0, .0}, {70, .0, .0}, {120, .0, .0}, 0.5, 5., 1.);
    EXPECT_EQ(particles[0], p1_expect);
    EXPECT_EQ(particles[1], p2_expect);
}

}  // namespace mol_sim
