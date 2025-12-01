#include "particles/boundaries/Outflow.h"

#include <gtest/gtest.h>

#include <limits>
#include <memory>
#include <numbers>

#include "io/outputWriter/XYZWriter.h"
#include "particles/Particle.h"
#include "particles/container/domain/Domain.h"
#include "physics/LennardJonesForce.h"
#include "testingUtils.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the VelocityReflection boundary condition
 * Base data/setup:
 * settings: default values below, other values set in tests
 * p:
 *  position: (5,5,5)
 *  velocity: (0,0,0)
 *  (mass: 1)
 * dimension: (10,10,10)
 * delta_t: 0.5
 * end_time = 0.5 (in most tests)
 */
class VelocityReflectionTest : public testing::Test {
   protected:
    SettingsParam settings;

    R3 p_x = {5.0, 5.0, 5.0};
    R3 p_v = {0.0, 0.0, 0.0};
    R3 dimension = {10.0, 10.0, 10.0};
    Particle p;
    double delta_t = 0.5;
    double end_time = 1 * delta_t;

    VelocityReflectionTest() : p{p_x, p_v, 1.0, 5., 1.} {}

    void SetUp() override {
        settings.domain = Domain(dimension);
        settings.delta_t = delta_t;
        settings.end_time = end_time;
        settings.cutoff = std::numeric_limits<double>::infinity();
    }
};
/**
 * @brief Tests that a particle is correctly reflected along the positive x boundary in a LinkedCellContainer.
 * */
TEST_F(VelocityReflectionTest, X_velocityreflect_linked) {
    p.getX() = {9.0, 5.0, 5.0};
    p.getV() = {4.0, 0.0, 4.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff);
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                               std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {9.0, 5.0, 7.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected along the negative x boundary in a SimpleContainer.
 * */
TEST_F(VelocityReflectionTest, X_velocityreflect_simple) {
    p.getX() = {1.0, 5.0, 5.0};
    p.getV() = {4.0, 0.0, 8.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {1.0, 5.0, 9.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected along the positive y boundary in a LinkedCellContainer.
 * */
TEST_F(VelocityReflectionTest, Y_velocityreflect_linked) {
    p.getX() = {5.0, 8.0, 5.0};
    p.getV() = {0.0, 4.0, 4.0};
    settings.end_time = 2 * delta_t;
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff);
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                               std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {5.0, 8.0, 9.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected along the negative y boundary in a SimpleContainer.
 * */
TEST_F(VelocityReflectionTest, Y_velocityreflect_simple) {
    p.getX() = {5.0, 1.0, 5.0};
    p.getV() = {4.0, 4.0, 0.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {7.0, 1.0, 5.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected along the positive z boundary in a LinkedCellContainer.
 * */
TEST_F(VelocityReflectionTest, Z_velocityreflect_linked) {
    p.getX() = {5.0, 5.0, 9.0};
    p.getV() = {0.0, -4.0, 4.0};
    settings.end_time = 2 * delta_t;
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff);
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                               std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {5.0, 1.0, 9.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected along the negative z boundary in a SimpleContainer.
 * */
TEST_F(VelocityReflectionTest, Z_velocityreflect_simple) {
    p.getX() = {5.0, 5.0, 1.0};
    p.getV() = {0.0, -2.0, -4.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {5.0, 4.0, 1.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected on a corner where two boundaries meet in a LinkedCellContainer.
 * */
TEST_F(VelocityReflectionTest, VelocityReflect_two_corner_linked) {
    p.getX() = {1.0, 1.0, 1.0};
    p.getV() = {-4.0, 4.0, 4.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff);
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {1.0, 2.0, 1.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected on a corner where two boundaries meet in a SimpleContainer.
 * */
TEST_F(VelocityReflectionTest, VelocityReflect_two_corner_simple) {
    p.getX() = {1.0, 1.0, 1.0};
    p.getV() = {-4.0, 4.0, 4.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {1.0, 2.0, 1.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected on a corner where three boundaries meet in a LinkedCellContainer.
 * */
TEST_F(VelocityReflectionTest, VelocityReflect_three_corner_linked) {
    p.getX() = {1.0, 9.0, 1.0};
    p.getV() = {-4.0, 4.0, 4.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff);
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {1.0, 9.0, 1.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected on a corner where three boundaries meet in a SimpleContainer.
 * */
TEST_F(VelocityReflectionTest, VelocityReflect_three_corner_simple) {
    p.getX() = {1.0, 9.0, 1.0};
    p.getV() = {-4.0, 4.0, 4.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {1.0, 9.0, 1.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected if shot diagonally in all directions in a LinkedCellContainer.
 * */
TEST_F(VelocityReflectionTest, VelocityReflect_diagonal_shot_linked) {
    p.getX() = {5.0, 9.0, 5.0};
    p.getV() = {4.0, 4.0, 4.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff);
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {6.0, 9.0, 6.0};
    EXPECT_EQ(p.getX(), expected);
}
/**
 * @brief Tests that a particle is correctly reflected if shot diagonally in all directions in a SimpleContainer.
 * */
TEST_F(VelocityReflectionTest, VelocityReflect_diagonal_shot_simple) {
    p.getX() = {5.0, 9.0, 5.0};
    p.getV() = {4.0, 4.0, 4.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    R3 expected = {6.0, 9.0, 6.0};
    EXPECT_EQ(p.getX(), expected);
}
}  // namespace mol_sim
