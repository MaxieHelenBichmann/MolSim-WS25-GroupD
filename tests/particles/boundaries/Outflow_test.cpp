#include <gtest/gtest.h>

#include <memory>
#include <numbers>
#include <limits>

#include "io/outputWriter/XYZWriter.h"
#include "particles/Particle.h"
#include "particles/boundaries/Outflow.h"
#include "particles/container/domain/Domain.h"
#include "physics/LennardJonesForce.h"
#include "testingUtils.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the Outflow boundary condition
 * Base data/setup:
 * settings: default values below, other values set in tests
 * p:
 *  position: (5,5,5)
 *  velocity: (0,0,0)
 *  (mass: 1)
 * dimension: (10,10,10)
 * delta_t: 0.5
 * end_time = 1
 */
class OutflowTest : public testing::Test {
   protected:
    SettingsParam settings;

    R3 p_x = {5.0, 5.0, 5.0};
    R3 p_v = {0.0, 0.0, 0.0};
    R3 dimension = {10.0, 10.0, 10.0};
    Particle p;
    double delta_t = 0.5;
    double end_time = 1 * delta_t;

    OutflowTest() : p{p_x, p_v, 1.0, 5., 1.} {}

    void SetUp() override {
        settings.domain = Domain(dimension);
        settings.delta_t = delta_t;
        settings.end_time = end_time;
        settings.cutoff = std::numeric_limits<double>::infinity();
        settings.setDefaults();
    }
};

/**
 * @brief Tests that outflowing particles in positive x are correctly erased in a LinkedCellContainer.
 * */
TEST_F(OutflowTest, X_outflow_linked) {
    p.getV() = {20.0, 0.0, 0.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                               std::make_unique<XYZWriter>());
    simulation.run();
    EXPECT_TRUE(particles.begin() == particles.end());
}
/**
 * @brief Tests that outflowing particles in negative x are correctly erased in a SimpleContainer.
 * */
TEST_F(OutflowTest, X_outflow_simple) {
    p.getV() = {-20.0, 0.0, 0.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    EXPECT_TRUE(particles.begin() == particles.end());
}
/**
 * @brief Tests that outflowing particles in positive y are correctly erased in a LinkedCellContainer.
 * */
TEST_F(OutflowTest, Y_outflow_linked) {
    p.getV() = {0.0, 20.0, 0.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                               std::make_unique<XYZWriter>());
    simulation.run();
    EXPECT_TRUE(particles.begin() == particles.end());
}
/**
 * @brief Tests that outflowing particles in negative y are correctly erased in a SimpleContainer.
 * */
TEST_F(OutflowTest, Y_outflow_simple) {
    p.getV() = {0.0, -20.0, 0.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    EXPECT_TRUE(particles.begin() == particles.end());
}
/**
 * @brief Tests that outflowing particles in positive z are correctly erased in a LinkedCellContainer.
 * */
TEST_F(OutflowTest, Z_outflow_linked) {
    p.getV() = {0.0, 0.0, 20.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                               std::make_unique<XYZWriter>());
    simulation.run();
    EXPECT_TRUE(particles.begin() == particles.end());
}
/**
 * @brief Tests that outflowing particles in negative z are correctly erased in a SimpleContainer.
 * */
TEST_F(OutflowTest, Z_outflow_simple) {
    p.getV() = {0.0, 0.0, -20.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    EXPECT_TRUE(particles.begin() == particles.end());
}
/**
 * @brief Tests that particles that are not OOB yet (but close to boundary) are NOT erased in a
 * LinkedCellContainer.
 * */
TEST_F(OutflowTest, outflow_linked_no_erase) {
    p.getV() = {9.0, 0.0, 0.0};
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, settings.cutoff.value());
    particles.addParticle(p);
    Simulation<LinkedCellContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                               std::make_unique<XYZWriter>());
    simulation.run();
    Particle expected({9.5, 5.0, 5.0}, p.getV(), p.getM(), p.getEpsilon(), p.getSigma(), p.getType());
    EXPECT_TRUE(particles[0] == expected);
}
/**
 * @brief Tests that particles that are not OOB yet (but close to boundary) are NOT erased in a
 * SimpleContainer.
 * */
TEST_F(OutflowTest, outflow_simple_no_erase) {
    p.getV() = {9.0, 0.0, 0.0};
    settings.container_type = "SIMPLE";
    SimpleContainer particles;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, std::make_unique<LennardJonesForce>(), settings,
                                           std::make_unique<XYZWriter>());
    simulation.run();
    Particle expected({9.5, 5.0, 5.0}, p.getV(), p.getM(), p.getEpsilon(), p.getSigma(), p.getType());
    EXPECT_TRUE(particles[0] == expected);
}
}  // namespace mol_sim
