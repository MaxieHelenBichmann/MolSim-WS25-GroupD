#include <gtest/gtest.h>

#include <numbers>
#include <memory>
#include <array>

#include "particles/Particle.h"
#include "utils/Settings.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Outflow.h"
#include "particles/container/domain/Domain.h"
#include "utils/Simulation.h"
#include "testingUtils.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the Reflecting boundary condition
 * Base data/setup:
 * settings: default values below, other values set in tests
 * p:
 *  position: (5,5,5)
 *  velocity: (0,0,0)
 *  (mass: 1)
 * dimension: (10,10,10)
 * delta_t: 0.1
 * end_time: 2 * delta_t = 0.2
 * boundaries: Reflecting, 
 */
class ReflectingTest : public testing::Test {
   protected:
    SettingsParam settings;
    
    R3 p_x = {5.0, 5.0, 5.0};
    R3 p_v = {0.0, 0.0, 0.0};
    R3 dimension = {10.0, 10.0, 10.0};
    Particle p;
    double delta_t = 0.1;
    double end_time = 2 * delta_t;

    ReflectingTest() : p{p_x, p_v, 1.0, 10., 6.} {}

    void SetUp() override {
        std::array<std::unique_ptr<Boundary>, 6> boundaries;
        boundaries[0] = std::make_unique<Reflecting>(BoundaryLocation::LEFT, false, 1, 20); 
        boundaries[1] = std::make_unique<Reflecting>(BoundaryLocation::RIGHT, false, 1, 20); 
        boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::UPPER, false, 1, 20); 
        boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::LOWER, false, 1, 20); 
        boundaries[4] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, false, 1, 20); 
        boundaries[5] = std::make_unique<Reflecting>(BoundaryLocation::BACK, false, 1, 20); 
        settings.domain = Domain(dimension, std::move(boundaries));
        settings.delta_t = delta_t;
        settings.end_time = end_time;
        settings.setDefaults();
    }
};
/**
 * @brief Tests that particles are correctly reflected on a x-boundary in a LinkedCellContainer. 
 * */
TEST_F(ReflectingTest, X_reflecting_linked) {
    p.getX() = {8.0, 5.0, 5.0};
    p.getV() = {15.0, 0.0, 0.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {6.2, 5.0, 5.0};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that particles are correctly reflected on a x-boundary in a SimpleContainer. 
 * */
TEST_F(ReflectingTest, X_reflecting_simple) {
    p.getX() = {2.0, 5.0, 5.0};
    p.getV() = {-15.0, 0.0, 0.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {3.8, 5.0, 5.0};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that particles are correctly reflected on a y-boundary in a LinkedCellContainer. 
 * */
TEST_F(ReflectingTest, Y_reflecting_linked) {
    p.getX() = {5.0, 8.0, 5.0};
    p.getV() = {0.0, 15.0, 0.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 6.2, 5.0};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that particles are correctly reflected on a y-boundary in a SimpleContainer. 
 * */
TEST_F(ReflectingTest, Y_reflecting_simple) {
    p.getX() = {5.0, 2.0, 5.0};
    p.getV() = {0.0, -15.0, 0.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 3.8, 5.0};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that particles are correctly reflected on a z-boundary in a LinkedCellContainer. 
 * */
TEST_F(ReflectingTest, Z_reflecting_linked) {
    p.getX() = {5.0, 5.0, 8.0};
    p.getV() = {0.0, 0.0, 15.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 5.0, 6.2};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that particles are correctly reflected on a z-boundary in a SimpleContainer. 
 * */
TEST_F(ReflectingTest, Z_reflecting_simple) {
    p.getX() = {5.0, 5.0, 2.0};
    p.getV() = {0.0, 0.0, -15.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 5.0, 3.8};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting particles at a corner of two boundaries works in a LinkedCellContainer.
 */
TEST_F(ReflectingTest, Reflecting_two_sided_corner_linked) {
    p.getX() = {2.0, 5.0, 2.0};
    p.getV() = {0.0, -15.0, -15.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {3.8, 5.0, 3.8};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting particles at a corner of three boundaries works in a LinkedCellContainer.
 */
TEST_F(ReflectingTest, Reflecting_three_sided_corner_linked) {
    p.getX() = {2.0, 2.0, 2.0};
    p.getV() = {-15.0, -15.0, -15.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {3.8, 3.8, 3.8};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting a particle that approaches a boundary at an angular trajectory 
 * (i.e. angle != 0 degrees) works in a LinkedCellContainer.
 */
TEST_F(ReflectingTest, Reflecting_Angular_trajectory_linked) {
    p.getX() = {2.0, 5.0, 2.0};
    p.getV() = {-15.0, 0.0, 15.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {3.8, 5.0, 5.0};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting a particle with a boundary taking in either the sigma or 
 * epsilon of the particle (or both) works in a LinkedCellContainer.
 */
TEST_F(ReflectingTest, Reflecting_particle_sigma_epsilon_linked) {
    p.getX() = {5.0, 5.0, 2.0};
    p.getV() = {0.0, 0.0, -15.0};
    p.setSigma(1.0);
    p.setEpsilon(20.0);
    std::array<std::unique_ptr<Boundary>, 6> boundaries;
    boundaries[0] = std::make_unique<Reflecting>(BoundaryLocation::LEFT, false, std::nullopt, std::nullopt); 
    boundaries[1] = std::make_unique<Reflecting>(BoundaryLocation::RIGHT, false, std::nullopt, std::nullopt); 
    boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::UPPER, false, std::nullopt, 20); 
    boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::LOWER, false, 1, std::nullopt); 
    boundaries[4] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, false, std::nullopt, std::nullopt); 
    boundaries[5] = std::make_unique<Reflecting>(BoundaryLocation::BACK, false, std::nullopt, std::nullopt); 
    settings.domain = Domain(dimension, std::move(boundaries));
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 5.0, 3.8};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting a particle while spawning ghost particles directly on the boundary 
 * works in a LinkedCellContainer.
 */
TEST_F(ReflectingTest, Reflecting_with_ghost_on_boundary_linked) {
    p.getX() = {5.0, 5.0, 2.5};
    p.getV() = {0.0, 0.0, -15.0};
    std::array<std::unique_ptr<Boundary>, 6> boundaries;
    boundaries[0] = std::make_unique<Reflecting>(BoundaryLocation::LEFT, true, 1, 20); 
    boundaries[1] = std::make_unique<Reflecting>(BoundaryLocation::RIGHT, true, 1, 20); 
    boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::UPPER, true, 1, 20); 
    boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::LOWER, true, 1, 20); 
    boundaries[4] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, true, 1, 20); 
    boundaries[5] = std::make_unique<Reflecting>(BoundaryLocation::BACK, true, 1, 20); 
    settings.domain = Domain(dimension, std::move(boundaries));
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<LinkedCellContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 5.0, 4.3};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting particles at a corner of two boundaries works in a SimpleContainer.
 */
TEST_F(ReflectingTest, Reflecting_two_sided_corner_simple) {
    p.getX() = {2.0, 5.0, 2.0};
    p.getV() = {0.0, -15.0, -15.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {3.8, 5.0, 3.8};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting particles at a corner of three boundaries works in a SimpleContainer.
 */
TEST_F(ReflectingTest, Reflecting_three_sided_corner_simple) {
    p.getX() = {2.0, 2.0, 2.0};
    p.getV() = {-15.0, -15.0, -15.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {3.8, 3.8, 3.8};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting a particle that approaches a boundary at an angular trajectory 
 * (i.e. angle != 0 degrees) works in a SimpleContainer.
 */
TEST_F(ReflectingTest, Reflecting_Angular_trajectory_simple) {
    p.getX() = {2.0, 5.0, 2.0};
    p.getV() = {-15.0, 0.0, 15.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {3.8, 5.0, 5.0};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting a particle with a boundary taking in either the sigma or 
 * epsilon of the particle (or both) works in a SimpleContainer.
 */
TEST_F(ReflectingTest, Reflecting_particle_sigma_epsilon_simple) {
    p.getX() = {5.0, 5.0, 2.0};
    p.getV() = {0.0, 0.0, -15.0};
    p.setSigma(1.0);
    p.setEpsilon(20.0);
    std::array<std::unique_ptr<Boundary>, 6> boundaries;
    boundaries[0] = std::make_unique<Reflecting>(BoundaryLocation::LEFT, false, std::nullopt, std::nullopt); 
    boundaries[1] = std::make_unique<Reflecting>(BoundaryLocation::RIGHT, false, std::nullopt, std::nullopt); 
    boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::UPPER, false, std::nullopt, 20); 
    boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::LOWER, false, 1, std::nullopt); 
    boundaries[4] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, false, std::nullopt, std::nullopt); 
    boundaries[5] = std::make_unique<Reflecting>(BoundaryLocation::BACK, false, std::nullopt, std::nullopt); 
    settings.domain = Domain(dimension, std::move(boundaries));
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 5.0, 3.8};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}
/**
 * @brief Tests that reflecting a particle while spawning ghost particles directly on the boundary 
 * works in a SimpleContainer.
 */
TEST_F(ReflectingTest, Reflecting_with_ghost_on_boundary_simple) {
    p.getX() = {5.0, 5.0, 2.5};
    p.getV() = {0.0, 0.0, -15.0};
    std::array<std::unique_ptr<Boundary>, 6> boundaries;
    boundaries[0] = std::make_unique<Reflecting>(BoundaryLocation::LEFT, true, 1, 20); 
    boundaries[1] = std::make_unique<Reflecting>(BoundaryLocation::RIGHT, true, 1, 20); 
    boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::UPPER, true, 1, 20); 
    boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::LOWER, true, 1, 20); 
    boundaries[4] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, true, 1, 20); 
    boundaries[5] = std::make_unique<Reflecting>(BoundaryLocation::BACK, true, 1, 20); 
    settings.domain = Domain(dimension, std::move(boundaries));
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    LennardJonesForce lj_force;
    Simulation<SimpleContainer, LennardJonesForce> simulation(particles, lj_force, settings);
    simulation.run();
    R3 expected = {5.0, 5.0, 4.3};
    EXPECT_TRUE(particles[0] == p);
    EXPECT_EQ(particles[0].getX(), expected);
}

}  // namespace mol_sim