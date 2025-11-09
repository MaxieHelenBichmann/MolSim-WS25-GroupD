#include "particles/container/SimpleContainer.h"

#include <gtest/gtest.h>

#include "particles/ParticleContainer.h"

namespace mol_sim {

static_assert(ParticleContainer<SimpleContainer>, "SimpleContainer must satisfy ParticleContainer concept");

/**
 * @brief Test Fixture for testing the SimpleContainer.
 *
 * Base Config/Data of the container is as follow:
 *
 * particles_empty = empty SimpleContainer
 * particles_full = SimpleContainer storing 4 Particles (content unimportant)
 *
 */
class SimpleContainerTest : public testing::Test {
   protected:
    SimpleContainer particles_empty;
    SimpleContainer particles_full;
    Particle p0;
    Particle p1;
    Particle p2;
    Particle p3;

    SimpleContainerTest()
        : particles_empty(), particles_full(), p0(Particle(0)), p1(Particle(1)), p2(Particle(2)), p3(Particle(3)) {}

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
TEST_F(SimpleContainerTest, testAccessSubscript) {
    EXPECT_TRUE(particles_full[0] == p0);
    EXPECT_TRUE(particles_full[2] == p2);
    EXPECT_TRUE(particles_full[3] == p3);
}

/**
 * @brief Tests correct behaviour of method size().
 */
TEST_F(SimpleContainerTest, testSize) {
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(0));
    particles_empty.addParticle(Particle(5));
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(4));
}

/**
 * @brief Tests correct behaviour of method empty().
 */
TEST_F(SimpleContainerTest, testEmpty) {
    EXPECT_FALSE(particles_full.empty());
    EXPECT_TRUE(particles_empty.empty());
    particles_empty.addParticle(Particle(4));
    EXPECT_FALSE(particles_empty.empty());
}

// modify

/**
 * @brief Tests correct behaviour of method clear().
 */
TEST_F(SimpleContainerTest, testClear) {
    ASSERT_EQ(particles_full.size(), static_cast<size_t>(4));
    particles_full.clear();
    EXPECT_TRUE(particles_full.empty());
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(0));
}

/**
 * @brief Tests correct behaviour of method reserve(size_t n).
 */
TEST_F(SimpleContainerTest, testReserve) {
    particles_empty.reserve(64);
    EXPECT_GE(particles_empty.capacity(), static_cast<size_t>(64));
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(0));
}

/**
 * @brief Tests correct behaviour of the method overload addParticle(Particle&& value).
 */
TEST_F(SimpleContainerTest, testAddParticleRval) {
    Particle tmp(6);
    particles_empty.addParticle(std::move(tmp));
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method overload addParticle(const Particle& value).
 */
TEST_F(SimpleContainerTest, testAddParticleConstLval) {
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
TEST_F(SimpleContainerTest, testAddParticleEmplaceNoType) {
    particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 5., 1.);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method overload
 * addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type).
 */
TEST_F(SimpleContainerTest, testAddParticleEmplace) {
    particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 8., 7., 69);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

// iterators

/**
 * @brief Tests correct behaviour of required non-const iterator begin().
 */
TEST_F(SimpleContainerTest, testBeginIterator) {
    auto full_it = particles_full.begin();
    auto full_it_end = particles_full.end();
    EXPECT_NE(full_it, full_it_end);
    EXPECT_EQ(&*full_it, particles_full.data());

    auto empty_it = particles_empty.begin();
    auto empty_it_end = particles_empty.end();
    EXPECT_EQ(empty_it, empty_it_end);
}

/**
 * @brief Tests correct behaviour of required non-const iterator end().
 */
TEST_F(SimpleContainerTest, testEndIterator) {
    auto n = static_cast<std::ptrdiff_t>(particles_full.size());
    EXPECT_EQ(&*(particles_full.end() - 1), particles_full.data() + (n - 1));  // NOLINT
}

/**
 * @brief Tests correct behaviour of required const iterators begin() and cbegin().
 */
TEST_F(SimpleContainerTest, testBeginConstInterator) {
    // begin()
    const SimpleContainer& cc = particles_full;
    auto cc_it = cc.begin();
    auto cc_it_end = cc.end();
    EXPECT_NE(cc_it, cc_it_end);
    EXPECT_EQ(&*cc_it, cc.data());

    // cbegin()
    auto it = particles_full.cbegin();
    auto it_end = particles_full.cend();
    EXPECT_NE(it, it_end);
    EXPECT_EQ(&*it, particles_full.data());
}

/**
 * @brief Tests correct behaviour of required const iterators end() and cend().
 */
TEST_F(SimpleContainerTest, testEndConstInterator) {
    // end()
    const SimpleContainer& cc = particles_full;
    auto n = static_cast<std::ptrdiff_t>(cc.size());
    EXPECT_EQ(&*(cc.end() - 1), cc.data() + (n - 1));  // NOLINT

    // cend()
    auto nc = static_cast<std::ptrdiff_t>(particles_full.size());
    EXPECT_EQ(&*(particles_full.cend() - 1), particles_full.data() + (nc - 1));  // NOLINT
}

// ParticleContainer: complex tests

/**
 * @brief Tests a sequence of read and write operations on a SimpleContainer.
 */
TEST(SimpleContainer, testSizeEmptyClearReserve) {
    SimpleContainer c;
    EXPECT_TRUE(c.empty());
    EXPECT_EQ(c.size(), static_cast<size_t>(0));

    c.reserve(64);
    EXPECT_GE(c.capacity(), static_cast<size_t>(64));

    c.addParticle(Particle(42));
    c.addParticle({1., 2., 3.}, {4., 5., 6.}, 1.0, 5., 1.);
    ASSERT_EQ(c.size(), static_cast<size_t>(2));

    c.clear();
    EXPECT_TRUE(c.empty());
    EXPECT_EQ(c.size(), static_cast<size_t>(0));
    EXPECT_EQ(c.begin(), c.end());
}

/**
 * @brief When we reserve a large amount of memory, we do not want to reallocate.
 */
TEST(SimpleContainer, testNoUnnecessaryReallocation) {
    SimpleContainer c;
    c.reserve(128);
    auto* old_data = c.data();
    for (int i = 0; i < 10; ++i) {
        c.addParticle(Particle(i));
    }
    EXPECT_EQ(old_data, c.data());
}

}  // namespace mol_sim
