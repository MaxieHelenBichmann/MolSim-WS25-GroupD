#include "utils/Simulation.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <io/outputWriter/XYZWriter.h>
#include <particles/container/ContainerRef.h>
#include <particles/container/domain/Domain.h>
#include <physics/pairwiseforces/GravitationalForce.h>
#include <physics/pairwiseforces/LennardJonesForce.h>

#include <limits>
#include <memory>

#include "io/CheckpointWriter.h"
#include "io/OutputWriter.h"
#include "io/StatsWriter.h"
#include "particles/Particle.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/SimpleContainer.h"
#include "testingUtils.h"
#include "utils/Vector.h"

namespace mol_sim {

// ==================== Mocks ====================

/**
 * @brief Mocks ForceSource.applyForce so we use better test values.
 */
class ForceMock : public PairwiseForceSource {
   public:
    // NOLINTNEXTLINE(bugprone-exception-escape)
    MOCK_METHOD(R3, applyForce, (const Particle& p1, const Particle& p2), (const, noexcept, override));
};

/**
 * @brief Mocks OutputWriter for testing without file I/O.
 */
class OutputWriterMock : public OutputWriter {
   public:
    MOCK_METHOD(void, plotParticles, (ContainerRef particles, const std::string& filename, int iteration),
                (const, override));
};

/**
 * @brief Mocks CheckpointWriter for testing without file I/O.
 *
 */
class CheckpointWriterMock : public CheckpointWriter {
   public:
    MOCK_METHOD(void, createCheckpoint,
                (SettingsParam & settings, const Domain& domain, ContainerRef particles, int iteration, size_t N),
                (const, override));
};

/**
 * @brief Mocks StatsWriter for testing without file I/O.
 *
 */
class StatsWriterMock : public StatsWriter {
   public:
    MOCK_METHOD(void, plotDiffusion, (ContainerRef particles, int iteration), (const));
    MOCK_METHOD(void, plotRDF, (ContainerRef particles, int iteration, const std::vector<Particle>& mirrored), (const));
};

// ==================== Container Traits ====================

template <typename Container>
struct ContainerTraits;

template <>
struct ContainerTraits<SimpleContainer> {
    static SimpleContainer create() { return SimpleContainer{}; }
};

template <>
struct ContainerTraits<LinkedCellContainer> {
    static LinkedCellContainer create() { return {{100.0, 100.0, 100.0}, std::numeric_limits<double>::infinity()}; }
};

// ==================== Base Test Fixture ====================

/**
 * @brief Typed test fixture for container-agnostic simulation tests.
 * @tparam Container The particle container type (SimpleContainer or LinkedCellContainer)
 */
template <typename Container>
class SimulationTestBase : public testing::Test {
   protected:
    Container particles = ContainerTraits<Container>::create();
    SettingsParam settings;
    // NOLINTNEXTLINE
    const double precision = 5e-8;

    void SetUp() override {
        settings.cutoff = std::numeric_limits<double>::infinity();
        settings.thermostat_freq = std::numeric_limits<size_t>::max();
        settings.thermo = false;
        settings.target_force_enabled = false;
        particles.clear();
    }
};

using ContainerTypes = testing::Types<SimpleContainer, LinkedCellContainer>;
TYPED_TEST_SUITE(SimulationTestBase, ContainerTypes);

template <typename Container>
class CalculateXTest : public SimulationTestBase<Container> {};
TYPED_TEST_SUITE(CalculateXTest, ContainerTypes);

template <typename Container>
class CalculateVTest : public SimulationTestBase<Container> {};
TYPED_TEST_SUITE(CalculateVTest, ContainerTypes);

template <typename Container>
class CalculateFTest : public SimulationTestBase<Container> {};
TYPED_TEST_SUITE(CalculateFTest, ContainerTypes);

template <typename Container>
class CalculateThermostatTest : public SimulationTestBase<Container> {};
TYPED_TEST_SUITE(CalculateThermostatTest, ContainerTypes);

template <typename Container>
class SimulationRunTest : public SimulationTestBase<Container> {};
TYPED_TEST_SUITE(SimulationRunTest, ContainerTypes);

/**
 * @brief Tests basic functionality of calculateX with respect to a
 delta_t != 1 and an initial velocity.
 */
TYPED_TEST(CalculateXTest, calculateX_dt) {
    this->settings.delta_t = 0.5;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 0.0};
    this->particles.addParticle(x, v, 1.0, 5., 1.);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<GravitationalForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateX();
    R3 expected = {5.0, -10.0, 0.0};
    EXPECT_EQ(this->particles[0].getX(), expected);
}

/**
 * @brief Tests calculateX when a particle has an initial velocity
 and an initial force that's acting on it.
 */
TYPED_TEST(CalculateXTest, calculateX_force) {
    this->settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, 0.0, 0.0};
    R3 f = {0.0, 10.0, 0.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    this->particles.addParticle(p);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<GravitationalForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateX();
    R3 expected = {10.0, 5.0, 0.0};
    EXPECT_EQ(this->particles[0].getX(), expected);
}

/**
 * @brief Tests basic functionality of calculateV for a single particle
 when no force is acting on it.
 */
TYPED_TEST(CalculateVTest, calculateV_no_force) {
    this->settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    this->particles.addParticle(x, v, 1.0, 5., 1.);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<GravitationalForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateV(1.);
    EXPECT_EQ(this->particles[0].getV(), v);
}

/**
 * @brief Tests basic functionalty of calculateV for a single particle
 when there is a force acting on it.
 */
TYPED_TEST(CalculateVTest, calculateV_simple) {
    this->settings.delta_t = 1;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 f = {1.0, 1.0, 1.0};
    Particle p(x, v, 1.0, 5., 1.);
    p.getF() = f;
    this->particles.addParticle(p);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<GravitationalForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateV(1.);
    R3 expected = {10.5, -19.5, 30.5};
    EXPECT_EQ(this->particles[0].getV(), expected);
}

/**
 * @brief Tests a more complex setting of calculateV where the
 particle has both a non-zero 'old_f' and 'f'.
 */
TYPED_TEST(CalculateVTest, calculateV_complex) {
    this->settings.delta_t = 0.5;
    R3 x = {0.0, 0.0, 0.0};
    R3 v = {10.0, -20.0, 30.0};
    R3 old_f = {-12.0, 15.0, -20.5};
    R3 f = {20.5, -10.0, 100.0};
    Particle p(x, v, 2.0, 5., 1.);
    p.getF() = f;
    p.getOldF() = old_f;
    this->particles.addParticle(p);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<GravitationalForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateV(1.);
    R3 expected = {11.0625, -19.375, 39.9375};
    EXPECT_EQ(this->particles[0].getV(), expected);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctlyf
 especially in regards to Newton's third law.
 The inter-particle force is 1 dimensional.
 */
TYPED_TEST(CalculateFTest, calculateF_simple2_pairwise) {
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.push_back(std::move(mock));
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateF(0);
    EXPECT_EQ(this->particles[0].getF(), f12);
    EXPECT_EQ(this->particles[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 2 particles are calculated correctly
 especially in regards to Newton's third law.
 The inter-particle force is 3 dimensional.
 */
TYPED_TEST(CalculateFTest, calculateF_complex2_pairwise) {
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    const Particle p2({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    R3 f12 = {102.52, -51.3, 135.711};
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.push_back(std::move(mock));
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateF(0);
    EXPECT_EQ(this->particles[0].getF(), f12);
    EXPECT_EQ(this->particles[1].getF(), -1.0 * f12);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly.
 Each inter-particle force is 1 dimensional.
 */
TYPED_TEST(CalculateFTest, calculateF_simple3_pairwise) {
    Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p2({2.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p3({3.0, 1.0, 1.0}, {.0, .0, .0}, {0., 0., 0.}, 1.0, 5., 1.);
    Particle p12({1.0, 1.0, 1.0}, {.0, .0, .0}, {10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p22({2.0, 1.0, 1.0}, {.0, .0, .0}, {-10.0, 0.0, 0.0}, 1.0, 5., 1.);
    Particle p32({3.0, 1.0, 1.0}, {.0, .0, .0}, {-20.0, 0.0, 0.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 0.0, 0.0};
    R3 f13 = {20.0, 0.0, 0.0};
    R3 f23 = {-10.0, 0.0, 0.0};
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    this->particles.addParticle(p3);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.push_back(std::move(mock));
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateF(0);
    R3 expected1 = {30.0, 0.0, 0.0};
    R3 expected2 = {-20.0, 0.0, 0.0};
    R3 expected3 = {-10.0, 0.0, 0.0};
    EXPECT_EQ(this->particles[0].getF(), expected1);
    EXPECT_EQ(this->particles[1].getF(), expected2);
    EXPECT_EQ(this->particles[2].getF(), expected3);
}

/**
 * @brief Tests that the forces of 3 particles are calculated correctly.
 Each of the inter-particle forces are 3 dimensional.
 */
TYPED_TEST(CalculateFTest, calculateF_complex3_pairwise) {
    const Particle p1({1.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p2({2.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    const Particle p3({3.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p12({1.0, 1.0, 1.0}, {.0, .0, .0}, {10.0, 5.0, 6.0}, 1.0, 5., 1.);
    Particle p22({2.0, 1.0, 1.0}, {.0, .0, .0}, {-10.0, -5.0, -6.0}, 1.0, 5., 1.);
    Particle p32({3.0, 1.0, 1.0}, {.0, .0, .0}, {-20.0, -10.0, -7.0}, 1.0, 5., 1.);
    R3 f12 = {10.0, 5.0, 6.0};
    R3 f13 = {20.0, 10.0, 7.0};
    R3 f23 = {-10.0, -5.0, -4.0};
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    this->particles.addParticle(p3);
    auto mock = std::make_unique<ForceMock>();
    EXPECT_CALL(*mock, applyForce(p1, p2)).Times(1).WillOnce(testing::Return(f12));
    EXPECT_CALL(*mock, applyForce(p12, p3)).Times(1).WillOnce(testing::Return(f13));
    EXPECT_CALL(*mock, applyForce(p22, p32)).Times(1).WillOnce(testing::Return(f23));
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.push_back(std::move(mock));
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.calculateF(0);
    R3 expected1 = {30.0, 15.0, 13.0};
    R3 expected2 = {-20.0, -10.0, -10.0};
    R3 expected3 = {-10.0, -5.0, -3.0};
    EXPECT_EQ(this->particles[0].getF(), expected1);
    EXPECT_EQ(this->particles[1].getF(), expected2);
    EXPECT_EQ(this->particles[2].getF(), expected3);
}

/**
 * @brief Tests that the Thermostat behaves correctly when the total energy of the system is 0.
 * There should not be a division by zero
 */
TYPED_TEST(CalculateThermostatTest, test_total_energy_0) {
    const Particle p{{1., 1., 0.}, {0., 0., 0.}, {0., 0., 0.}, 1., 5., 1.};
    this->particles.addParticle(p);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    XYZWriter writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.getTotalEnergy() = 0;
    EXPECT_EQ(simulation.calculateThermostatFactor(), 1.);
}

/**
 * @brief Tests the beta is computed correctly when there is no delta T limit applied.
 */
TYPED_TEST(CalculateThermostatTest, test_factor_without_delta) {
    const Particle p1{{1., 1., 0.}, {10., 0., 0.}, {0., 0., 0.}, 1., 5., 1.};
    const Particle p2{{2., 1., 0.}, {5., 5., 0.}, {0., 0., 0.}, 3., 5., 1.};
    const Particle p3{{3., 1., 0.}, {3., 4., 5.}, {0., 0., 0.}, 2., 5., 1.};
    // total energy should be 175
    // current temp should be 350/9 = 38.88...
    double expected = 1.133893419;
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    this->particles.addParticle(p3);
    this->settings.target_temp = 50;
    this->settings.delta_temp = std::numeric_limits<double>::infinity();
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    XYZWriter writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.getTotalEnergy() = 175;
    EXPECT_NEAR(simulation.calculateThermostatFactor(), expected, this->precision);
}

/**
 * @brief Tests that the thermostat behaves correctly when the system is already at the target temperature.
 * I.E. no change in velocities should happen.
 */
TYPED_TEST(CalculateThermostatTest, test_factor_at_target_temp) {
    const Particle p1{{1., 1., 0.}, {10., 0., 0.}, {0., 0., 0.}, 1., 5., 1.};
    const Particle p2{{2., 1., 0.}, {5., 5., 0.}, {0., 0., 0.}, 3., 5., 1.};
    // total energy at 125
    // current temp at 125/3
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    this->settings.target_temp = 125. / 3.;
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    XYZWriter writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.getTotalEnergy() = 125;
    EXPECT_NEAR(simulation.calculateThermostatFactor(), 1., this->precision);
}

/**
 * @brief Tests that the limit in temperature change specified by delta temp is applied correctly.
 */
TYPED_TEST(CalculateThermostatTest, test_delta_temp) {
    const Particle p1{{1., 1., 0.}, {10., 0., 0.}, {0., 0., 0.}, 1., 5., 1.};
    const Particle p2{{2., 1., 0.}, {5., 5., 0.}, {0., 0., 0.}, 3., 5., 1.};
    const Particle p3{{3., 1., 0.}, {3., 4., 5.}, {0., 0., 0.}, 2., 5., 1.};
    // total energy should be 175
    // current temp should be 350/9 = 38.88...
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    this->particles.addParticle(p3);
    this->settings.target_temp = 50;
    this->settings.delta_temp = 1.;
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    XYZWriter writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.getTotalEnergy() = 175;
    double b = simulation.calculateThermostatFactor();

    double new_energy = 0.;

    for (Particle& p : this->particles) {
        p.getV() = p.getV() * b;
        new_energy += p.getM() * R3::scalarProduct(p.getV(), p.getV());
    }
    new_energy *= 0.5;

    double new_temp = (2. * new_energy) / (this->settings.dimensions * this->particles.size());
    double curr_temp = 350. / 9;
    EXPECT_NEAR(new_temp - curr_temp, 1., this->precision);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 gravitational forces and 2 particles.
 */
TYPED_TEST(SimulationRunTest, run_gravitational_timestep) {
    this->settings.delta_t = 0.5;
    this->settings.end_time = 0.5;
    this->settings.target_temp = 100. / 3.;
    Particle p1({2.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({1.0, 1.0, 1.0}, {10.0, .0, .0}, {.0, .0, .0}, 0.5, 5., 1.);
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<GravitationalForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.run();
    // After 1 timestep: positions update first (with old_f=0), then forces calculated, then velocities
    R3 p1_x_expect = {2.0, 1.0, 1.0};
    R3 p2_x_expect = {6.0, 1.0, 1.0};
    R3 p1_v_expect = {0.0078125, 0.0, 0.0};
    R3 p2_v_expect = {9.984375, 0.0, 0.0};
    EXPECT_R3_NEAR(this->particles[0].getX(), p1_x_expect, this->precision);
    EXPECT_R3_NEAR(this->particles[1].getX(), p2_x_expect, this->precision);
    EXPECT_R3_NEAR(this->particles[0].getV(), p1_v_expect, this->precision);
    EXPECT_R3_NEAR(this->particles[1].getV(), p2_v_expect, this->precision);
}

/**
 * @brief Tests that a single timestep is calculated correctly in run() with
 Lennard-Jones forces and 2 particles.
 */
TYPED_TEST(SimulationRunTest, run_lennardjones_timestep) {
    this->settings.delta_t = 0.5;
    this->settings.end_time = 0.5;
    this->settings.target_temp = 0.;
    // Place particles at distance 1 apart for strong LJ interaction
    Particle p1({5.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    Particle p2({4.0, 1.0, 1.0}, {.0, .0, .0}, {.0, .0, .0}, 1.0, 5., 1.);
    this->particles.addParticle(p1);
    this->particles.addParticle(p2);
    std::vector<std::unique_ptr<PairwiseForceSource>> pairwise_forces;
    pairwise_forces.emplace_back(std::make_unique<LennardJonesForce>());
    std::vector<std::unique_ptr<SingleForceSource>> single_forces;
    OutputWriterMock writer;
    CheckpointWriterMock checkpoint_writer;
    StatsWriterMock stat_writer;
    Simulation<TypeParam> simulation(this->particles, pairwise_forces, single_forces, this->settings, writer,
                                     checkpoint_writer, stat_writer);
    simulation.run();
    // At distance=1 with sigma=1, epsilon=5: F = -24*5*(1-2)*direction = 120*direction
    // p1-p2 = (1,0,0), so F on p1 = (120,0,0), F on p2 = (-120,0,0)
    // x updates with old_f=0, so x stays same
    // Then new forces calculated at same positions: F = (120,0,0) on p1
    // v1 = 0 + 0.5*0.5/1*(0+120,0,0) = (30,0,0)
    // v2 = 0 + 0.5*0.5/1*(0-120,0,0) = (-30,0,0)
    R3 p1_x_expect = {5.0, 1.0, 1.0};
    R3 p2_x_expect = {4.0, 1.0, 1.0};
    R3 p1_v_expect = {30.0, 0.0, 0.0};
    R3 p2_v_expect = {-30.0, 0.0, 0.0};

    EXPECT_R3_NEAR(this->particles[0].getX(), p1_x_expect, this->precision);
    EXPECT_R3_NEAR(this->particles[1].getX(), p2_x_expect, this->precision);
    EXPECT_R3_NEAR(this->particles[0].getV(), p1_v_expect, this->precision);
    EXPECT_R3_NEAR(this->particles[1].getV(), p2_v_expect, this->precision);
}

}  // namespace mol_sim
