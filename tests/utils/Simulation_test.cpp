#include "utils/Simulation.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <utils/Settings.h>

#include <utility>

#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

class ForceMock : public ForceSource {
   public:
    MOCK_METHOD(R3, applyForce, (const Particle& p1, const Particle& p2), (const, override));
};

class OutputWriterMock : public OutputWriter {
   public:
    MOCK_METHOD(void, plotParticles, (ContainerRef particles, const std::string& filename, int iteration), (override));
};

class SimulationTest : public testing::Test {
   protected:
    SimpleContainer particles;
    SettingsParam settings;

   public:
    void SetUp() override {
        settings.setDefaults();
        particles.clear();
    }
};

TEST_F(SimulationTest, calculateX_dt) {
    settings.delta_t = 0.5;

    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 0.0};
    particles.addParticle(x, v, 1.0);
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.calculateX();
    R3 expected = {5.0, -10.0, 0.0};
    EXPECT_EQ(particles[0].getX(), expected);
}

TEST_F(SimulationTest, calculateX_force) {
    settings.delta_t = 1;
    settings.force_type = GRAVITATIONAL;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, 0.0, 0.0};
    R3 f = {0.0, 10.0, 0.0};
    Particle p(x, v, 1.0);
    p.getF() = f;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.calculateX();
    R3 expected = {10.0, 5.0, 0.0};
    EXPECT_EQ(particles[0].getX(), expected);
}

TEST_F(SimulationTest, calculateV_no_force) {
    settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    particles.addParticle(x, v, 1.0);
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.calculateV();
    EXPECT_EQ(particles[0].getV(), v);
}

TEST_F(SimulationTest, calculateV_simple) {
    settings.delta_t = 1;
    settings.force_type = GRAVITATIONAL;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 f = {1.0, 1.0, 1.0};
    Particle p(x, v, 1.0);
    p.getF() = f;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.calculateV();
    R3 expected = {10.5, -19.5, 30.5};
    EXPECT_EQ(particles[0].getV(), expected);
}

TEST_F(SimulationTest, calculateV_complex) {
    settings.delta_t = 0.5;
    settings.force_type = GRAVITATIONAL;
    settings.setDefaults();
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 old_f = {-12.0, 15.0, -20.5};
    R3 f = {20.5, -10.0, 100.0};
    Particle p(x, v, 2.0);
    p.getF() = f;
    p.getOldF() = old_f;
    particles.addParticle(p);
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.calculateV();
    R3 expected = {11.0625, -19.375, 39.9375};
    EXPECT_EQ(particles[0].getV(), expected);
}

TEST_F(SimulationTest, calculateF_simple2_pairwise) {
    auto force = std::make_unique<ForceMock>();
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({.0, .0, .0}, {.0, .0, .0}, 1.0);
    const Particle p2({.0, .0, .0}, {.0, .0, .0}, 1.0);
    R3 f12 = {10.0, 0.0, 0.0};
    EXPECT_CALL(*force, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    particles.addParticle(p1);
    particles.addParticle(p2);
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.setForceSource(std::move(force));
    simulation.calculateF();
    EXPECT_EQ(particles[0].getF(), f12);
    EXPECT_EQ(particles[1].getF(), -1.0 * f12);
}

TEST_F(SimulationTest, calculateF_complex2_pairwise) {
    auto force = std::make_unique<ForceMock>();
    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({.0, .0, .0}, {.0, .0, .0}, 1.0);
    const Particle p2({.0, .0, .0}, {.0, .0, .0}, 1.0);
    R3 f12 = {102.52, -51.3, 135.711};
    EXPECT_CALL(*force, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    particles.addParticle(p1);
    particles.addParticle(p2);
    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.setForceSource(std::move(force));
    simulation.calculateF();
    EXPECT_EQ(particles[0].getF(), f12);
    EXPECT_EQ(particles[1].getF(), -1.0 * f12);
}

TEST_F(SimulationTest, calculateF_simple3_pairwise) {
    auto force = std::make_unique<ForceMock>();
    // initial positions and velocities irrelevant since we're mocking
    Particle p1({1.0, .0, .0}, {.0, .0, .0}, 1.0);
    Particle p2({2.0, .0, .0}, {.0, .0, .0}, 1.0);
    Particle p3({3.0, .0, .0}, {.0, .0, .0}, 1.0);

    Particle p12({1.0, .0, .0}, {.0, .0, .0}, 1.0);
    p12.getF() = {10.0, 0.0, 0.0};

    Particle p22({2.0, .0, .0}, {.0, .0, .0}, 1.0);
    p22.getF() = {-10.0, 0.0, 0.0};

    Particle p32({3.0, .0, .0}, {.0, .0, .0}, 1.0);
    p32.getF() = {-20.0, 0.0, 0.0};

    R3 f12 = {10.0, 0.0, 0.0};
    R3 f13 = {20.0, 0.0, 0.0};
    R3 f23 = {-10.0, 0.0, 0.0};

    EXPECT_CALL(*force, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*force, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*force, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));

    particles.addParticle(p1);
    particles.addParticle(p2);
    particles.addParticle(p3);

    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.setForceSource(std::move(force));
    simulation.calculateF();

    R3 expected1 = {30.0, 0.0, 0.0};
    R3 expected2 = {-20.0, 0.0, 0.0};
    R3 expected3 = {-10.0, 0.0, 0.0};

    EXPECT_EQ(particles[0].getF(), expected1);
    EXPECT_EQ(particles[1].getF(), expected2);
    EXPECT_EQ(particles[2].getF(), expected3);
}

TEST_F(SimulationTest, calculateF_complex3_pairwise) {
    auto force = std::make_unique<ForceMock>();

    // initial positions and velocities irrelevant since we're mocking
    const Particle p1({1.0, .0, .0}, {.0, .0, .0}, 1.0);
    const Particle p2({2.0, .0, .0}, {.0, .0, .0}, 1.0);
    const Particle p3({3.0, .0, .0}, {.0, .0, .0}, 1.0);

    Particle p12({1.0, .0, .0}, {.0, .0, .0}, 1.0);
    p12.getF() = {10.0, 5.0, 6.0};

    Particle p22({2.0, .0, .0}, {.0, .0, .0}, 1.0);
    p22.getF() = {-10.0, -5.0, -6.0};

    Particle p32({3.0, .0, .0}, {.0, .0, .0}, 1.0);
    p32.getF() = {-20.0, -10.0, -7.0};

    R3 f12 = {10.0, 5.0, 6.0};
    R3 f13 = {20.0, 10.0, 7.0};
    R3 f23 = {-10.0, -5.0, -4.0};

    EXPECT_CALL(*force, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*force, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*force, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));

    particles.addParticle(p1);
    particles.addParticle(p2);
    particles.addParticle(p3);

    Simulation<SimpleContainer> simulation(particles, settings);
    simulation.setForceSource(std::move(force));
    simulation.calculateF();

    R3 expected1 = {30.0, 15.0, 13.0};
    R3 expected2 = {-20.0, -10.0, -10.0};
    R3 expected3 = {-10.0, -5.0, -3.0};

    EXPECT_EQ(particles[0].getF(), expected1);
    EXPECT_EQ(particles[1].getF(), expected2);
    EXPECT_EQ(particles[2].getF(), expected3);
}

// TODO: test run(). 2 tests, one for simple GRAVITATIONAL (solar system or sum), the other for simple LENNARDJONES,
// each 3 iterations or something, nothing crazy.
}  // namespace mol_sim
