#include "utils/Simulation.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <particles/container/domain/Domain.h>
#include <physics/GravitationalForce.h>
#include <physics/LennardJonesForce.h>

#include <limits>

#include "io/OutputWriter.h"
#include "particles/Particle.h"
#include "particles/container/SimpleContainer.h"
#include "particles/container/LinkedCellContainer.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Mocks ForceSource.applyForce so we use better test values.
 *
 */
class ForceMock : public ForceSource {
   public:
    // Mocked method for applying force between two particles.
    // NOLINTNEXTLINE(bugprone-exception-escape)
    MOCK_METHOD(R3, applyForce, (const Particle& p1, const Particle& p2), (const, noexcept, override));
};

/**
 * @brief Mocks OutputWriter for testing without file I/O.
 *
 */
class OutputWriterMock : public OutputWriter {
   public:
    MOCK_METHOD(void, plotParticles, (ContainerRef particles, const std::string& filename, int iteration),
                (const, override));
};

/**
 * @brief Test fixture for the following tests. Contains two attributes
 that will be used throughout the tests.
 1) 'particles_simple' which is a SimpleContainer containing the particles we perform the SimpleContainer tests on
 2) 'particles_linked' which is a LinkedCellContainer containing the particles we perform the LinkedCellContainer tests on
 3) 'domain_size' which represents the domain_size (important for particles_linked)
 4) 'settings' which is a SettingsParam for setting simulation parameters
 */
class SimulationTest : public testing::Test {
   protected:
    SimpleContainer particles_simple;
    R3 domain_size = {100.0, 100.0, 100.0};
    LinkedCellContainer* particles_linked_heap_ptr = new LinkedCellContainer(domain_size, std::numeric_limits<double>::infinity());
    LinkedCellContainer particles_linked = *particles_linked_heap_ptr;
    SettingsParam settings;

   public:
    void SetUp() override {
        settings.cutoff = std::numeric_limits<double>::infinity();
        particles_simple.clear();
        particles_linked.clear();
    }

    void TearDown() override {
        delete particles_linked_heap_ptr;
    }
};

/**
 * @brief Tests basic functionality of calculateX with respect to a
 delta_t != 1 and an initial velocity in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, calculateX_dt_simplecontainer) {
    settings.delta_t = 0.5;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 0.0};
    particles_simple.addParticle(x, v, 1.0, 5., 1.);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *force_source, settings, *writer);
    simulation.calculateX();
    R3 expected = {5.0, -10.0, 0.0};
    EXPECT_EQ(particles_simple[0].getX(), expected);
}

/**
 * @brief Tests calculateX when a particle has an initial velocity
 and an initial force that's acting on it in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, calculateX_force_simplecontainer) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, 0.0, 0.0};
    R3 f = {0.0, 10.0, 0.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    particles_simple.addParticle(p);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *force_source, settings, *writer);
    simulation.calculateX();
    R3 expected = {10.0, 5.0, 0.0};
    EXPECT_EQ(particles_simple[0].getX(), expected);
}

/**
 * @brief Tests basic functionality of calculateV for a single particle
 when no force is acting on it in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, calculateV_no_force_simplecontainer) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    particles_simple.addParticle(x, v, 1.0, 5., 1.);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *force_source, settings, *writer);
    simulation.calculateV();
    EXPECT_EQ(particles_simple[0].getV(), v);
}

/**
 * @brief Tests basic functionalty of calculateV for a single particle
 when there is a force acting on it in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, calculateV_simple_simplecontainer) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 f = {1.0, 1.0, 1.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    particles_simple.addParticle(p);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *force_source, settings, *writer);
    simulation.calculateV();
    R3 expected = {10.5, -19.5, 30.5};
    EXPECT_EQ(particles_simple[0].getV(), expected);
}

/**
 * @brief Tests a more complex setting of calculateV where the
 particle has both a non-zero 'old_f' and 'f' in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, calculateV_complex_simplecontainer) {
    settings.delta_t = 0.5;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 old_f = {-12.0, 15.0, -20.5};
    R3 f = {20.5, -10.0, 100.0};
    Particle p(x, v, 2.0, 5., 1.);
    p.getF() = f;
    p.getOldF() = old_f;
    particles_simple.addParticle(p);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *force_source, settings, *writer);
    simulation.calculateV();
    R3 expected = {11.0625, -19.375, 39.9375};
    EXPECT_EQ(particles_simple[0].getV(), expected);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctly
 especially in regards to Newton's third law in a SimpleContainer. 
 The inter-particle force is 1 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_simple2_pairwise_simplecontainer) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    particles_simple.addParticle(p1);
    particles_simple.addParticle(p2);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *mock, settings, *writer);
    simulation.calculateF();
    EXPECT_EQ(particles_simple[0].getF(), f12);
    EXPECT_EQ(particles_simple[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctly
 especially in regards to Newton's third law in a SimpleContainer. 
 The inter-*mock is 3 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_complex2_pairwise_simplecontainer) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {102.52, -51.3, 135.711};
    particles_simple.addParticle(p1);
    particles_simple.addParticle(p2);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *mock, settings, *writer);
    simulation.calculateF();
    EXPECT_EQ(particles_simple[0].getF(), f12);
    EXPECT_EQ(particles_simple[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly.
 Each inter-particle force is 1 dimensional in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, calculateF_simple3_pairwise_simplecontainer) {
    // initial positions and velocities irrelevant since we're mocking
    Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p2({2.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p3({3.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p12({1.0, 1.0, 1.0}, {.0, .0, .0}, {10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p22({2.0, 1.0, 1.0}, {.0, .0, .0}, {-10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p32({3.0, 1.0, 1.0}, {.0, .0, .0}, {-20.0, 0.0, 0.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    R3 f13 = {20.0, 0.0, 0.0};
    R3 f23 = {-10.0, 0.0, 0.0};
    particles_simple.addParticle(p1);
    particles_simple.addParticle(p2);
    particles_simple.addParticle(p3);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *mock, settings, *writer);
    simulation.calculateF();
    R3 expected1 = {30.0, 0.0, 0.0};
    R3 expected2 = {-20.0, 0.0, 0.0};
    R3 expected3 = {-10.0, 0.0, 0.0};
    EXPECT_EQ(particles_simple[0].getF(), expected1);
    EXPECT_EQ(particles_simple[1].getF(), expected2);
    EXPECT_EQ(particles_simple[2].getF(), expected3);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly.
 Each of the inter-particle forces are 3 dimensional in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, calculateF_complex3_pairwise_simplecontainer) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p2({2.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p3({3.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p12({1.0, 1.0, 1.0}, {.0, .0, .0}, {10.0, 5.0, 6.0}, 1.0, 5., 1.);
    Particle p22({2.0, 1.0, 1.0}, {.0, .0, .0}, {-10.0, -5.0, -6.0}, 1.0, 5., 1.);
    Particle p32({3.0, 1.0, 1.0}, {.0, .0, .0}, {-20.0, -10.0, -7.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 5.0, 6.0};
    R3 f13 = {20.0, 10.0, 7.0};
    R3 f23 = {-10.0, -5.0, -4.0};
    particles_simple.addParticle(p1);
    particles_simple.addParticle(p2);
    particles_simple.addParticle(p3);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *mock, settings, *writer);

    simulation.calculateF();
    R3 expected1 = {30.0, 15.0, 13.0};
    R3 expected2 = {-20.0, -10.0, -10.0};
    R3 expected3 = {-10.0, -5.0, -3.0};
    EXPECT_EQ(particles_simple[0].getF(), expected1);
    EXPECT_EQ(particles_simple[1].getF(), expected2);
    EXPECT_EQ(particles_simple[2].getF(), expected3);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 gravitational forces and 2 particles in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, run_gravitational_timestep_simplecontainer) {
    settings.delta_t = 0.5;
    settings.end_time = 0.5;
    Particle p1({2.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({1.0, 1.0, 1.0}, {10.0, .0, .0}, {.0, .0, .0}, 0.5, 5., 1.);
    particles_simple.addParticle(p1);
    particles_simple.addParticle(p2);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *force_source, settings, *writer);
    simulation.run();
    // After 1 timestep: positions update first (with old_f=0), then forces calculated, then velocities
    Particle p1_expect({2.0, 1.0, 1.0}, {0.0078125, .0, .0}, {0.03125, .0, .0}, 1.0, 5., 1.);
    Particle p2_expect({6.0, 1.0, 1.0}, {9.984375, .0, .0}, {-0.03125, .0, .0}, 0.5, 5., 1.);
    p1_expect.getOldX() = {2.0, 1.0, 1.0};
    p2_expect.getOldX() = {1.0, 1.0, 1.0};
    EXPECT_EQ(particles_simple[0], p1_expect);
    EXPECT_EQ(particles_simple[1], p2_expect);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 Lennard-Jones forces and 2 particles in a SimpleContainer.
 *
 */
TEST_F(SimulationTest, run_lennardjones_timestep_simplecontainer) {
    settings.delta_t = 0.5;
    settings.end_time = 0.5;
    // Place particles_simple at distance 1 apart for strong LJ interaction
    Particle p1({5.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({4.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    particles_simple.addParticle(p1);
    particles_simple.addParticle(p2);
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<SimpleContainer> simulation(particles_simple, *force_source, settings, *writer);
    simulation.run();
    // At distance=1 with sigma=1, epsilon=5: F = -24*5*(1-2)*direction = 120*direction
    // p1-p2 = (1,0,0), so F on p1 = (120,0,0), F on p2 = (-120,0,0)
    // x updates with old_f=0, so x stays same
    // Then new forces calculated at same positions: F = (120,0,0) on p1
    // v1 = 0 + 0.5*0.5/1*(0+120,0,0) = (30,0,0)
    // v2 = 0 + 0.5*0.5/1*(0-120,0,0) = (-30,0,0)
    Particle p1_expect({5.0, 1.0, 1.0}, {30.0, .0, .0}, {120.0, .0, .0}, 1.0, 5., 1.);
    Particle p2_expect({4.0, 1.0, 1.0}, {-30.0, .0, .0}, {-120.0, .0, .0}, 1.0, 5., 1.);
    p1_expect.getOldX() = {5.0, 1.0, 1.0};
    p2_expect.getOldX() = {4.0, 1.0, 1.0};
    EXPECT_EQ(particles_simple[0], p1_expect);
    EXPECT_EQ(particles_simple[1], p2_expect);
}

/**
 * @brief Tests basic functionality of calculateX with respect to a
 delta_t != 1 and an initial velocity in a LinkedCellContainer.
 *
 */
TEST_F(SimulationTest, calculateX_dt_linkedcellcontainer) {
    settings.delta_t = 0.5;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 0.0};
    particles_linked.addParticle(x, v, 1.0, 5., 1.);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *force_source, settings, *writer);
    simulation.calculateX();
    R3 expected = {5.0, -10.0, 0.0};
    EXPECT_EQ(particles_linked[0].getX(), expected);
}

/**
 * @brief Tests calculateX when a particle has an initial velocity
 and an initial force that's acting on it in a LinkedCellContainer.
 *
 */
TEST_F(SimulationTest, calculateX_force_linkedcellcontainer) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, 0.0, 0.0};
    R3 f = {0.0, 10.0, 0.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    particles_linked.addParticle(p);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *force_source, settings, *writer);
    simulation.calculateX();
    R3 expected = {10.0, 5.0, 0.0};
    EXPECT_EQ(particles_linked[0].getX(), expected);
}

/**
 * @brief Tests basic functionality of calculateV for a single particle
 when no force is acting on it in a LinkedCellContainer.
 *
 */
TEST_F(SimulationTest, calculateV_no_force_linkedcellcontainer) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    particles_linked.addParticle(x, v, 1.0, 5., 1.);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *force_source, settings, *writer);
    simulation.calculateV();
    EXPECT_EQ(particles_linked[0].getV(), v);
}

/**
 * @brief Tests basic functionalty of calculateV for a single particle
 when there is a force acting on it in a LinkedCellContainer.
 *
 */
TEST_F(SimulationTest, calculateV_simple_linkedcellcontainer) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 f = {1.0, 1.0, 1.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    particles_linked.addParticle(p);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *force_source, settings, *writer);
    simulation.calculateV();
    R3 expected = {10.5, -19.5, 30.5};
    EXPECT_EQ(particles_linked[0].getV(), expected);
}

/**
 * @brief Tests a more complex setting of calculateV where the
 particle has both a non-zero 'old_f' and 'f' in a LinkedCellContainer.
 *
 */
TEST_F(SimulationTest, calculateV_complex_linkedcellcontainer) {
    settings.delta_t = 0.5;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 old_f = {-12.0, 15.0, -20.5};
    R3 f = {20.5, -10.0, 100.0};
    Particle p(x, v, 2.0, 5., 1.);
    p.getF() = f;
    p.getOldF() = old_f;
    particles_linked.addParticle(p);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *force_source, settings, *writer);
    simulation.calculateV();
    R3 expected = {11.0625, -19.375, 39.9375};
    EXPECT_EQ(particles_linked[0].getV(), expected);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctly
 especially in regards to Newton's third law in a LinkedCellContainer. 
 The inter-particle force is 1 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_simple2_linkedcellcontainer) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    particles_linked.addParticle(p1);
    particles_linked.addParticle(p2);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *mock, settings, *writer);
    simulation.calculateF();
    EXPECT_EQ(particles_linked[0].getF(), f12);
    EXPECT_EQ(particles_linked[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctly
 especially in regards to Newton's third law in a LinkedCellContainer. 
 The inter-mock is 3 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_complex2_linkedcellcontainer) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {102.52, -51.3, 135.711};
    particles_linked.addParticle(p1);
    particles_linked.addParticle(p2);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *mock, settings, *writer);
    simulation.calculateF();
    EXPECT_EQ(particles_linked[0].getF(), f12);
    EXPECT_EQ(particles_linked[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly in a LinkedCellContainer.
 Each inter-particle force is 1 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_simple3_linkedcellcontainer) {
    // initial positions and velocities irrelevant since we're mocking
    Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p2({2.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p3({3.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p12({1.0, 1.0, 1.0}, {.0, .0, .0}, {10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p22({2.0, 1.0, 1.0}, {.0, .0, .0}, {-10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p32({3.0, 1.0, 1.0}, {.0, .0, .0}, {-20.0, 0.0, 0.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    R3 f13 = {20.0, 0.0, 0.0};
    R3 f23 = {-10.0, 0.0, 0.0};
    particles_linked.addParticle(p1);
    particles_linked.addParticle(p2);
    particles_linked.addParticle(p3);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *mock, settings, *writer);
    simulation.calculateF();
    R3 expected1 = {30.0, 0.0, 0.0};
    R3 expected2 = {-20.0, 0.0, 0.0};
    R3 expected3 = {-10.0, 0.0, 0.0};
    EXPECT_EQ(particles_linked[0].getF(), expected1);
    EXPECT_EQ(particles_linked[1].getF(), expected2);
    EXPECT_EQ(particles_linked[2].getF(), expected3);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly in a LinkedCellContainer.
 Each of the inter-particle forces are 3 dimensional.
 *
 */
TEST_F(SimulationTest, calculateF_complex3_pairwise_linkedcellcontainer) {
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p2({2.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p3({3.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p12({1.0, 1.0, 1.0}, {.0, .0, .0}, {10.0, 5.0, 6.0}, 1.0, 5., 1.);
    Particle p22({2.0, 1.0, 1.0}, {.0, .0, .0}, {-10.0, -5.0, -6.0}, 1.0, 5., 1.);
    Particle p32({3.0, 1.0, 1.0}, {.0, .0, .0}, {-20.0, -10.0, -7.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 5.0, 6.0};
    R3 f13 = {20.0, 10.0, 7.0};
    R3 f23 = {-10.0, -5.0, -4.0};
    particles_linked.addParticle(p1);
    particles_linked.addParticle(p2);
    particles_linked.addParticle(p3);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *mock, settings, *writer);
    simulation.calculateF();
    R3 expected1 = {30.0, 15.0, 13.0};
    R3 expected2 = {-20.0, -10.0, -10.0};
    R3 expected3 = {-10.0, -5.0, -3.0};
    EXPECT_EQ(particles_linked[0].getF(), expected1);
    EXPECT_EQ(particles_linked[1].getF(), expected2);
    EXPECT_EQ(particles_linked[2].getF(), expected3);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 gravitational forces and 2 particles in a LinkedCellContainer.
 *
 */
TEST_F(SimulationTest, run_gravitational_timestep_linkedcellcontainer) {
    settings.delta_t = 0.5;
    settings.end_time = 0.5;
    Particle p1({2.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({1.0, 1.0, 1.0}, {10.0, .0, .0}, {.0, .0, .0}, 0.5, 5., 1.);
    particles_linked.addParticle(p1);
    particles_linked.addParticle(p2);
    auto force_source = std::make_unique<GravitationalForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *force_source, settings, *writer);
    simulation.run();
    // After 1 timestep: positions update first (with old_f=0), then forces calculated, then velocities
    Particle p1_expect({2.0, 1.0, 1.0}, {0.0078125, .0, .0}, {0.03125, .0, .0}, 1.0, 5., 1.);
    Particle p2_expect({6.0, 1.0, 1.0}, {9.984375, .0, .0}, {-0.03125, .0, .0}, 0.5, 5., 1.);
    p1_expect.getOldX() = {2.0, 1.0, 1.0};
    p2_expect.getOldX() = {1.0, 1.0, 1.0};
    EXPECT_EQ(particles_linked[0], p1_expect);
    EXPECT_EQ(particles_linked[1], p2_expect);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 Lennard-Jones forces and 2 particles in a LinkedCellContainer.
 *
 */
TEST_F(SimulationTest, run_lennardjones_timestep_linkedcellcontainer) {
    settings.delta_t = 0.5;
    settings.end_time = 0.5;
    // Place particles at distance 1 apart for strong LJ interaction
    Particle p1({5.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({4.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    particles_linked.addParticle(p1);
    particles_linked.addParticle(p2);
    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<OutputWriterMock>();
    Simulation<LinkedCellContainer> simulation(particles_linked, *force_source, settings, *writer);
    simulation.run();
    // At distance=1 with sigma=1, epsilon=5: F = -24*5*(1-2)*direction = 120*direction
    // p1-p2 = (1,0,0), so F on p1 = (120,0,0), F on p2 = (-120,0,0)
    // x updates with old_f=0, so x stays same
    // Then new forces calculated at same positions: F = (120,0,0) on p1
    // v1 = 0 + 0.5*0.5/1*(0+120,0,0) = (30,0,0)
    // v2 = 0 + 0.5*0.5/1*(0-120,0,0) = (-30,0,0)
    Particle p1_expect({5.0, 1.0, 1.0}, {30.0, .0, .0}, {120.0, .0, .0}, 1.0, 5., 1.);
    Particle p2_expect({4.0, 1.0, 1.0}, {-30.0, .0, .0}, {-120.0, .0, .0}, 1.0, 5., 1.);
    p1_expect.getOldX() = {5.0, 1.0, 1.0};
    p2_expect.getOldX() = {4.0, 1.0, 1.0};
    EXPECT_EQ(particles_linked[0], p1_expect);
    EXPECT_EQ(particles_linked[1], p2_expect);
}

}  // namespace mol_sim
