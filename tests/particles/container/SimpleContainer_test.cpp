#include "particles/container/SimpleContainer.h"

#include <gtest/gtest.h>

#include "particles/ParticleContainer.h"

namespace mol_sim {

// Note: Basic ParticleContainer interface tests are in ParticleContainer_test.cpp
// This file contains SimpleContainer-specific tests only.

static_assert(ParticleContainer<SimpleContainer>, "SimpleContainer must satisfy ParticleContainer concept");

/**
 * @brief Test Fixture for SimpleContainer-specific tests.
 */
class SimpleContainerSpecificTest : public testing::Test {
   protected:
    SimpleContainer particles_empty;
    SimpleContainer particles_full;
    Particle p0;
    Particle p1;
    Particle p2;
    Particle p3;

    SimpleContainerSpecificTest()
        : particles_empty({10.0, 10.0, 10.0}, 2.5),
          particles_full({10.0, 10.0, 10.0}, 2.5),
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

// ============================================================================
// SimpleContainer-specific tests (vector-based behavior)
// ============================================================================

/**
 * @brief Tests that reserve increases capacity.
 */
TEST_F(SimpleContainerSpecificTest, ReserveIncreasesCapacity) {
    particles_empty.reserve(64);
    EXPECT_GE(particles_empty.capacity(), 64U);
    EXPECT_EQ(particles_empty.size(), 0U);
}

/**
 * @brief Tests that subscript access returns particles in order.
 */
TEST_F(SimpleContainerSpecificTest, SubscriptOrderPreserved) {
    EXPECT_TRUE(particles_full[0] == p0);
    EXPECT_TRUE(particles_full[1] == p1);
    EXPECT_TRUE(particles_full[2] == p2);
    EXPECT_TRUE(particles_full[3] == p3);
}

/**
 * @brief Tests that erase maintains order for remaining elements.
 */
TEST_F(SimpleContainerSpecificTest, ErasePreservesOrder) {
    auto it = particles_full.begin();
    ++it;  // points to p1
    particles_full.eraseParticle(it);
    EXPECT_EQ(particles_full.size(), 3U);
    EXPECT_TRUE(particles_full[0] == p0);
    EXPECT_TRUE(particles_full[1] == p2);
    EXPECT_TRUE(particles_full[2] == p3);
}

/**
 * @brief Tests data() method returns pointer to first element.
 */
TEST_F(SimpleContainerSpecificTest, DataPointer) {
    auto full_it = particles_full.begin();
    EXPECT_EQ(&*full_it, particles_full.data());
}

/**
 * @brief Tests end iterator arithmetic.
 */
TEST_F(SimpleContainerSpecificTest, EndIteratorArithmetic) {
    auto n = static_cast<std::ptrdiff_t>(particles_full.size());
    EXPECT_EQ(&*(particles_full.end() - 1), particles_full.data() + (n - 1));  // NOLINT
}

/**
 * @brief Tests const data pointer.
 */
TEST_F(SimpleContainerSpecificTest, ConstDataPointer) {
    const SimpleContainer& cc = particles_full;
    auto cc_it = cc.begin();
    EXPECT_EQ(&*cc_it, cc.data());

    auto it = particles_full.cbegin();
    EXPECT_EQ(&*it, particles_full.data());
}

/**
 * @brief Tests const end iterator arithmetic.
 */
TEST_F(SimpleContainerSpecificTest, ConstEndIteratorArithmetic) {
    const SimpleContainer& cc = particles_full;
    auto n = static_cast<std::ptrdiff_t>(cc.size());
    EXPECT_EQ(&*(cc.end() - 1), cc.data() + (n - 1));  // NOLINT

    auto nc = static_cast<std::ptrdiff_t>(particles_full.size());
    EXPECT_EQ(&*(particles_full.cend() - 1), particles_full.data() + (nc - 1));  // NOLINT
}

/**
 * @brief When we reserve a large amount of memory, we do not want to reallocate.
 */
TEST(SimpleContainerSpecific, NoUnnecessaryReallocation) {
    SimpleContainer c({10.0, 10.0, 10.0}, 2.5);
    c.reserve(128);
    auto* old_data = c.data();
    for (int i = 0; i < 10; ++i) {
        c.addParticle(Particle(i));
    }
    EXPECT_EQ(old_data, c.data());
}

/**
 * @brief Tests proximity iterator with infinite radius returns all particles.
 */
TEST(SimpleContainerSpecific, ProximityIteratorInfiniteRadius) {
    SimpleContainer particles_inf({10.0, 10.0, 10.0}, std::numeric_limits<double>::infinity());
    R3 v{0.0, 0.0, 0.0};
    particles_inf.addParticle(R3{2.0, 2.0, 2.0}, v, 1.0, 1.0, 1.0);
    particles_inf.addParticle(R3{4.0, 4.0, 4.0}, v, 1.0, 1.0, 1.0);
    particles_inf.addParticle(R3{6.0, 6.0, 6.0}, v, 1.0, 1.0, 1.0);
    R3 center{3.0, 3.0, 3.0};

    auto it = particles_inf.proximityBegin(center, particles_inf.size());
    auto end = particles_inf.proximityEnd(center);
    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_EQ(count, 3U);
}

}  // namespace mol_sim
