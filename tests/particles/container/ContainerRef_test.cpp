#include "particles/container/ContainerRef.h"

#include <gtest/gtest.h>

#include "particles/ParticleContainer.h"
#include "particles/container/SimpleContainer.h"

namespace mol_sim {

static_assert(ParticleContainer<ContainerRef>, "ContainerRef must satisfy ParticleContainer concept");

/**
 * @brief Test Fixture for testing the ContainerRef.
 *
 * Base Config/Data of the container is as follow:
 *
 * particles_empty = empty ContainerRef with a SimpleContainer as storage
 * particles_full = ContainerRef storing 4 Particles (content unimportant) with a SimpleContainer as storage
 *
 */
class ContainerRefTest : public testing::Test {
   protected:
    SimpleContainer storage1;
    SimpleContainer storage2;
    ContainerRef particles_empty;
    ContainerRef particles_full;
    Particle p0;
    Particle p1;
    Particle p2;
    Particle p3;

    ContainerRefTest()
        : storage1(),
          storage2(),
          particles_empty(storage1),
          particles_full(storage2),
          p0(Particle(0)),
          p1(Particle(1)),
          p2(Particle(2)),
          p3(Particle(3)) {}

    void SetUp() override {
        particles_empty.clear();
        particles_full.clear();
        particles_full.addParticle(p0);
        particles_full.addParticle(p1);
        particles_full.addParticle(p2);
        particles_full.addParticle(p3);
    }
};

// retrieve data

/**
 * @brief Tests correct access with the subscript operator.
 */
TEST_F(ContainerRefTest, testAccessSubscript) {
    EXPECT_TRUE(particles_full[0] == p0 || particles_full[0] == p1 || particles_full[0] == p2 ||
                particles_full[0] == p3);
    EXPECT_TRUE(particles_full[2] == p0 || particles_full[2] == p1 || particles_full[2] == p2 ||
                particles_full[2] == p3);
}

/**
 * @brief Tests correct behaviour of method size().
 */
TEST_F(ContainerRefTest, testSize) {
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(0));
    particles_empty.addParticle(Particle(5));
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(4));
}

/**
 * @brief Tests correct behaviour of method empty().
 */
TEST_F(ContainerRefTest, testEmpty) {
    EXPECT_FALSE(particles_full.empty());
    EXPECT_TRUE(particles_empty.empty());
    particles_empty.addParticle(Particle(4));
    EXPECT_FALSE(particles_empty.empty());
}

// modify

/**
 * @brief Tests correct behaviour of method clear().
 */
TEST_F(ContainerRefTest, testClear) {
    ASSERT_EQ(particles_full.size(), static_cast<size_t>(4));
    particles_full.clear();
    EXPECT_TRUE(particles_full.empty());
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(0));
}

/**
 * @brief Tests correct behaviour of method reserve(size_t n).
 */
TEST_F(ContainerRefTest, testReserve) {
    particles_empty.reserve(64);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(0));
}

/**
 * @brief Tests correct behaviour of the method overload addParticle(Particle&& value).
 */
TEST_F(ContainerRefTest, testAddParticleRval) {
    Particle tmp(6);
    particles_empty.addParticle(std::move(tmp));
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method overload addParticle(const Particle& value).
 */
TEST_F(ContainerRefTest, testAddParticleConstLval) {
    const Particle tmp(7);
    particles_empty.addParticle(tmp);
    particles_full.addParticle(tmp);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(5));
    EXPECT_TRUE(particles_empty[0] == tmp);
    EXPECT_TRUE(particles_full[4] == tmp);
}

/**
 * @brief Tests correct behaviour of the method overload addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg,
 * double m_arg).
 */
TEST_F(ContainerRefTest, testAddParticleEmplaceNoType) {
    particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 5., 1.);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method overload
 * addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type).
 */
TEST_F(ContainerRefTest, testAddParticleEmplace) {
    particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 8., 9., 69);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method updateParticlePosition.
 */
TEST_F(ContainerRefTest, testUpdateParticlePosition) {
    R3 new_x = {9.0, 8.0, 7.0};
    auto it = particles_full.begin();
    particles_full.updateParticlePosition(it, new_x);
    EXPECT_EQ(particles_full.begin()->getX(), new_x);
}

// iterators

/**
 * @brief Tests correct behaviour of required non-const iterator begin().
 */
TEST_F(ContainerRefTest, testBeginIterator) {
    EXPECT_NE(particles_full.begin(), particles_full.end());
    EXPECT_EQ(particles_empty.begin(), particles_empty.end());
}

/**
 * @brief Tests correct behaviour of required const iterators begin() and cbegin().
 */
TEST_F(ContainerRefTest, testBeginConstInterator) {
    // begin()
    const ContainerRef& cc = particles_full;
    EXPECT_NE(cc.begin(), cc.end());

    // cbegin()
    EXPECT_NE(particles_full.cbegin(), particles_full.cend());
}

// ParticleContainer: complex tests

/**
 * @brief Tests a sequence of read and write operations on a ContainerRef.
 */
TEST(ContainerRef, testSizeEmptyClearReserve) {
    SimpleContainer s;
    ContainerRef c(s);
    EXPECT_TRUE(c.empty());
    EXPECT_EQ(c.size(), static_cast<size_t>(0));

    c.reserve(64);
    EXPECT_EQ(c.size(), static_cast<size_t>(0));

    c.addParticle(Particle(42));
    c.addParticle({1., 2., 3.}, {4., 5., 6.}, 1.0, 5., 1.);
    ASSERT_EQ(c.size(), static_cast<size_t>(2));

    c.clear();
    EXPECT_TRUE(c.empty());
    EXPECT_EQ(c.size(), static_cast<size_t>(0));
    EXPECT_EQ(c.begin(), c.end());
}

}  // namespace mol_sim
