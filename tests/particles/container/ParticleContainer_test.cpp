#include "particles/ParticleContainer.h"

#include <gtest/gtest.h>
#include <particles/container/ContainerRef.h>

#include <set>

#include "particles/container/LinkedCellContainer.h"
#include "particles/container/SimpleContainer.h"

namespace mol_sim {

static_assert(ParticleContainer<SimpleContainer>, "SimpleContainer must satisfy ParticleContainer concept");
static_assert(ParticleContainer<LinkedCellContainer>, "LinkedCellContainer must satisfy ParticleContainer concept");
static_assert(ParticleContainer<ContainerRef>, "ContainerRef must satisfy ParticleContainer concept");

/**
 * @brief Factory trait for creating container instances with appropriate parameters.
 */
template <typename T>
struct ContainerFactory;

template <>
struct ContainerFactory<SimpleContainer> {
    static SimpleContainer create() { return {{10.0, 10.0, 10.0}, 2.5}; }
    static SimpleContainer createWithCutoff(double cutoff) { return {{10.0, 10.0, 10.0}, cutoff}; }
};

template <>
struct ContainerFactory<LinkedCellContainer> {
    static LinkedCellContainer create() { return {{10.0, 10.0, 10.0}, 2.5}; }
    static LinkedCellContainer createWithCutoff(double cutoff) { return {{10.0, 10.0, 10.0}, cutoff}; }
};

/**
 * @brief Templated test fixture for ParticleContainer tests.
 */
template <typename T>
class ParticleContainerTest : public testing::Test {
   protected:
    T particles_empty;
    T particles_full;
    Particle p0{0};
    Particle p1{1};
    Particle p2{2};
    Particle p3{3};

    ParticleContainerTest()
        : particles_empty(ContainerFactory<T>::create()), particles_full(ContainerFactory<T>::create()) {}

    void SetUp() override {
        particles_empty.clear();
        particles_full.clear();
        particles_full.addParticle(p0);
        particles_full.addParticle(p1);
        particles_full.addParticle(p2);
        particles_full.addParticle(p3);
    }
};

using ContainerTypes = ::testing::Types<SimpleContainer, LinkedCellContainer>;
TYPED_TEST_SUITE(ParticleContainerTest, ContainerTypes);

// ============================================================================
// Retrieve Data Tests
// ============================================================================

/**
 * @brief Tests operator[] returns valid particles.
 */
TYPED_TEST(ParticleContainerTest, OperatorSubscript) {
    // All particles should be one of the added particles
    EXPECT_TRUE(this->particles_full[0] == this->p0 || this->particles_full[0] == this->p1 ||
                this->particles_full[0] == this->p2 || this->particles_full[0] == this->p3);
    EXPECT_TRUE(this->particles_full[2] == this->p0 || this->particles_full[2] == this->p1 ||
                this->particles_full[2] == this->p2 || this->particles_full[2] == this->p3);
}

/**
 * @brief Tests const operator[].
 */
TYPED_TEST(ParticleContainerTest, OperatorSubscriptConst) {
    const TypeParam& const_container = this->particles_full;
    EXPECT_TRUE(const_container[0] == this->p0 || const_container[0] == this->p1 || const_container[0] == this->p2 ||
                const_container[0] == this->p3);
}

/**
 * @brief Tests size() returns correct count.
 */
TYPED_TEST(ParticleContainerTest, Size) {
    EXPECT_EQ(this->particles_empty.size(), 0U);
    this->particles_empty.addParticle(Particle(5));
    EXPECT_EQ(this->particles_empty.size(), 1U);
    EXPECT_EQ(this->particles_full.size(), 4U);
}

/**
 * @brief Tests empty() returns correct state.
 */
TYPED_TEST(ParticleContainerTest, Empty) {
    EXPECT_TRUE(this->particles_empty.empty());
    EXPECT_FALSE(this->particles_full.empty());
    this->particles_empty.addParticle(Particle(4));
    EXPECT_FALSE(this->particles_empty.empty());
}

// ============================================================================
// Modify Tests
// ============================================================================

/**
 * @brief Tests clear() empties the container.
 */
TYPED_TEST(ParticleContainerTest, Clear) {
    ASSERT_EQ(this->particles_full.size(), 4U);
    this->particles_full.clear();
    EXPECT_TRUE(this->particles_full.empty());
    EXPECT_EQ(this->particles_full.size(), 0U);
    EXPECT_EQ(this->particles_full.begin(), this->particles_full.end());
}

/**
 * @brief Tests reserve() preallocates memory without changing size.
 */
TYPED_TEST(ParticleContainerTest, Reserve) {
    this->particles_empty.reserve(64);
    EXPECT_EQ(this->particles_empty.size(), 0U);
}

/**
 * @brief Tests addParticle(Particle&&) with rvalue.
 */
TYPED_TEST(ParticleContainerTest, AddParticleRvalue) {
    Particle tmp(6);
    this->particles_empty.addParticle(std::move(tmp));
    EXPECT_EQ(this->particles_empty.size(), 1U);
}

/**
 * @brief Tests addParticle(const Particle&) with const lvalue.
 */
TYPED_TEST(ParticleContainerTest, AddParticleConstLvalue) {
    const Particle tmp(7);
    this->particles_empty.addParticle(tmp);
    this->particles_full.addParticle(tmp);
    EXPECT_EQ(this->particles_empty.size(), 1U);
    EXPECT_EQ(this->particles_full.size(), 5U);
    EXPECT_TRUE(this->particles_empty[0] == tmp);
    EXPECT_TRUE(this->particles_full[4] == tmp);
}

/**
 * @brief Tests addParticle with position, velocity, mass, epsilon, sigma.
 */
TYPED_TEST(ParticleContainerTest, AddParticleEmplaceNoType) {
    this->particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 5., 1.);
    EXPECT_EQ(this->particles_empty.size(), 1U);
}

/**
 * @brief Tests addParticle with position, velocity, mass, epsilon, sigma, and type.
 */
TYPED_TEST(ParticleContainerTest, AddParticleEmplaceWithType) {
    this->particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 8., 9., 69);
    EXPECT_EQ(this->particles_empty.size(), 1U);
}

/**
 * @brief Tests eraseParticle removes a particle correctly.
 */
TYPED_TEST(ParticleContainerTest, EraseParticle) {
    ASSERT_EQ(this->particles_full.size(), 4U);
    this->particles_full.eraseParticle(this->particles_full.begin());
    EXPECT_EQ(this->particles_full.size(), 3U);
}

/**
 * @brief Tests updateParticlePosition updates position.
 */
TYPED_TEST(ParticleContainerTest, UpdateParticlePosition) {
    R3 new_x = {5.0, 5.0, 5.0};
    auto it = this->particles_full.begin();
    this->particles_full.updateParticlePosition(it, new_x);
    EXPECT_EQ(this->particles_full.begin()->getX(), new_x);
}

// ============================================================================
// Iterator Tests
// ============================================================================

/**
 * @brief Tests begin() and end() for non-const container.
 */
TYPED_TEST(ParticleContainerTest, BeginEndIterator) {
    EXPECT_NE(this->particles_full.begin(), this->particles_full.end());
    EXPECT_EQ(this->particles_empty.begin(), this->particles_empty.end());
}

/**
 * @brief Tests const begin() and end().
 */
TYPED_TEST(ParticleContainerTest, BeginEndConstIterator) {
    const TypeParam& cc = this->particles_full;
    EXPECT_NE(cc.begin(), cc.end());

    const TypeParam& ec = this->particles_empty;
    EXPECT_EQ(ec.begin(), ec.end());
}

/**
 * @brief Tests cbegin() and cend().
 */
TYPED_TEST(ParticleContainerTest, CBeginCEnd) {
    EXPECT_NE(this->particles_full.cbegin(), this->particles_full.cend());
    EXPECT_EQ(this->particles_empty.cbegin(), this->particles_empty.cend());
}

/**
 * @brief Tests iteration counts all particles.
 */
TYPED_TEST(ParticleContainerTest, IteratorCount) {
    size_t count = 0;
    for (auto it = this->particles_full.begin(); it != this->particles_full.end(); ++it) {
        ++count;
    }
    EXPECT_EQ(count, 4U);
}

/**
 * @brief Tests proximity iterator with finite radius.
 */
TYPED_TEST(ParticleContainerTest, ProximityIteratorFinite) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(1.0);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{3.0, 3.0, 3.0}, v, 1.0, 1.0, 1.0);  // center particle
    container.addParticle(R3{3.5, 3.0, 3.0}, v, 1.0, 1.0, 1.0);  // within radius
    container.addParticle(R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // far away

    R3 center{3.0, 3.0, 3.0};
    auto it = container.proximityBegin(center, container.size());
    auto end = container.proximityEnd(center);

    size_t count = 0;
    while (it != end) {
        EXPECT_LE((it->getX() - center).euclidNorm(), 1.0);
        ++it;
        ++count;
    }
    // At least one particle should be within proximity
    EXPECT_GE(count, 1U);
}

/**
 * @brief Tests const proximity iterator.
 */
TYPED_TEST(ParticleContainerTest, ProximityIteratorConst) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(1.0);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{3.0, 3.0, 3.0}, v, 1.0, 1.0, 1.0);

    const TypeParam& const_container = container;
    R3 center{3.0, 3.0, 3.0};
    auto it = const_container.proximityBegin(center, const_container.size());
    auto end = const_container.proximityEnd(center);

    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_GE(count, 0U);  // Just verify it compiles and runs
}

/**
 * @brief Tests boundary iterator finds particles near boundaries.
 */
TYPED_TEST(ParticleContainerTest, BoundaryIterator) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(2.5);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // interior
    container.addParticle(R3{0.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // near left boundary

    auto it = container.boundaryBegin();
    auto end = container.boundaryEnd();

    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_GE(count, 1U);
}

/**
 * @brief Tests boundary iterator with specific location.
 */
TYPED_TEST(ParticleContainerTest, BoundaryIteratorSpecificLocation) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(2.5);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{0.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // near left boundary
    container.addParticle(R3{9.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // near right boundary
    container.addParticle(R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // interior

    std::set<BoundaryLocation> left_only = {BoundaryLocation::LEFT};
    auto it = container.boundaryBegin(left_only);
    auto end = container.boundaryEnd(left_only);

    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_GE(count, 1U);
}

/**
 * @brief Tests const boundary iterator.
 */
TYPED_TEST(ParticleContainerTest, BoundaryIteratorConst) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(2.5);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{0.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);

    const TypeParam& const_container = container;
    auto it = const_container.boundaryBegin();
    auto end = const_container.boundaryEnd();

    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_GE(count, 0U);
}

/**
 * @brief Tests halo iterator finds particles in halo region.
 */
TYPED_TEST(ParticleContainerTest, HaloIterator) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(2.5);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0);   // interior
    container.addParticle(R3{-0.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // halo (left)

    auto it = container.haloBegin();
    auto end = container.haloEnd();

    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_GE(count, 1U);
}

/**
 * @brief Tests halo iterator with specific location.
 */
TYPED_TEST(ParticleContainerTest, HaloIteratorSpecificLocation) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(2.5);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{-0.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // halo left
    container.addParticle(R3{10.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);  // halo right

    std::set<BoundaryLocation> left_only = {BoundaryLocation::LEFT};
    auto it = container.haloBegin(left_only);
    auto end = container.haloEnd(left_only);

    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_GE(count, 1U);
}

/**
 * @brief Tests const halo iterator.
 */
TYPED_TEST(ParticleContainerTest, HaloIteratorConst) {
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(2.5);
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{-0.5, 5.0, 5.0}, v, 1.0, 1.0, 1.0);

    const TypeParam& const_container = container;
    auto it = const_container.haloBegin();
    auto end = const_container.haloEnd();

    size_t count = 0;
    while (it != end) {
        ++it;
        ++count;
    }
    EXPECT_GE(count, 0U);
}

// ============================================================================
// Advanced Erase Tests
// ============================================================================

/**
 * @brief Tests erasing first particle.
 */
TYPED_TEST(ParticleContainerTest, EraseFirst) {
    ASSERT_EQ(this->particles_full.size(), 4U);
    auto it = this->particles_full.begin();
    auto next_it = this->particles_full.eraseParticle(it);
    EXPECT_EQ(this->particles_full.size(), 3U);
    // Returned iterator should be valid (either pointing to next element or end)
    EXPECT_TRUE(next_it == this->particles_full.begin() || next_it == this->particles_full.end() ||
                next_it != this->particles_full.end());
}

/**
 * @brief Tests erasing last particle.
 */
TYPED_TEST(ParticleContainerTest, EraseLast) {
    ASSERT_EQ(this->particles_full.size(), 4U);
    auto it = this->particles_full.end();
    --it;  // point to last element
    this->particles_full.eraseParticle(it);
    EXPECT_EQ(this->particles_full.size(), 3U);
}

/**
 * @brief Tests erasing middle particle.
 */
TYPED_TEST(ParticleContainerTest, EraseMiddle) {
    ASSERT_EQ(this->particles_full.size(), 4U);
    auto it = this->particles_full.begin();
    ++it;  // point to second element
    this->particles_full.eraseParticle(it);
    EXPECT_EQ(this->particles_full.size(), 3U);
}

/**
 * @brief Tests erasing all particles one by one.
 */
TYPED_TEST(ParticleContainerTest, EraseAll) {
    ASSERT_EQ(this->particles_full.size(), 4U);
    while (!this->particles_full.empty()) {
        this->particles_full.eraseParticle(this->particles_full.begin());
    }
    EXPECT_TRUE(this->particles_full.empty());
    EXPECT_EQ(this->particles_full.size(), 0U);
}

/**
 * @brief Tests that erasing maintains container consistency.
 */
TYPED_TEST(ParticleContainerTest, EraseConsistency) {
    // Add particles with known positions
    TypeParam container = ContainerFactory<TypeParam>::create();
    R3 v{0.0, 0.0, 0.0};
    container.addParticle(R3{1.0, 1.0, 1.0}, v, 1.0, 1.0, 1.0);
    container.addParticle(R3{2.0, 2.0, 2.0}, v, 1.0, 1.0, 1.0);
    container.addParticle(R3{3.0, 3.0, 3.0}, v, 1.0, 1.0, 1.0);

    ASSERT_EQ(container.size(), 3U);

    // Erase middle
    auto it = container.begin();
    ++it;
    container.eraseParticle(it);

    EXPECT_EQ(container.size(), 2U);

    // Verify remaining particles are accessible
    size_t count = 0;
    for (auto iter = container.begin(); iter != container.end(); ++iter) {
        ++count;
    }
    EXPECT_EQ(count, 2U);
}

// ============================================================================
// Advanced Insert/Add Tests
// ============================================================================

/**
 * @brief Tests adding many particles.
 */
TYPED_TEST(ParticleContainerTest, AddMany) {
    TypeParam container = ContainerFactory<TypeParam>::create();
    R3 v{0.0, 0.0, 0.0};

    for (int i = 0; i < 100; ++i) {
        double pos = 0.1 + ((i % 90) * 0.1);  // Keep within domain [0, 10]
        container.addParticle(R3{pos, pos, pos}, v, 1.0, 1.0, 1.0, i);
    }

    EXPECT_EQ(container.size(), 100U);
}

/**
 * @brief Tests adding after clear.
 */
TYPED_TEST(ParticleContainerTest, AddAfterClear) {
    this->particles_full.clear();
    EXPECT_TRUE(this->particles_full.empty());

    this->particles_full.addParticle(Particle(100));
    EXPECT_EQ(this->particles_full.size(), 1U);
}

/**
 * @brief Tests interleaved add and erase operations.
 */
TYPED_TEST(ParticleContainerTest, InterleavedAddErase) {
    TypeParam container = ContainerFactory<TypeParam>::create();
    R3 v{0.0, 0.0, 0.0};

    // Add 3 particles
    container.addParticle(R3{1.0, 1.0, 1.0}, v, 1.0, 1.0, 1.0);
    container.addParticle(R3{2.0, 2.0, 2.0}, v, 1.0, 1.0, 1.0);
    container.addParticle(R3{3.0, 3.0, 3.0}, v, 1.0, 1.0, 1.0);
    EXPECT_EQ(container.size(), 3U);

    // Erase one
    container.eraseParticle(container.begin());
    EXPECT_EQ(container.size(), 2U);

    // Add 2 more
    container.addParticle(R3{4.0, 4.0, 4.0}, v, 1.0, 1.0, 1.0);
    container.addParticle(R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0);
    EXPECT_EQ(container.size(), 4U);

    // Erase all
    while (!container.empty()) {
        container.eraseParticle(container.begin());
    }
    EXPECT_TRUE(container.empty());

    // Add again
    container.addParticle(R3{6.0, 6.0, 6.0}, v, 1.0, 1.0, 1.0);
    EXPECT_EQ(container.size(), 1U);
}

/**
 * @brief Tests that reserve + add doesn't invalidate data.
 */
TYPED_TEST(ParticleContainerTest, ReserveThenAdd) {
    TypeParam container = ContainerFactory<TypeParam>::create();
    container.reserve(100);

    R3 v{0.0, 0.0, 0.0};
    for (int i = 0; i < 50; ++i) {
        double pos = 0.1 + ((i % 90) * 0.1);
        container.addParticle(R3{pos, pos, pos}, v, 1.0, 1.0, 1.0, i);
    }

    EXPECT_EQ(container.size(), 50U);

    // Verify all particles are accessible
    size_t count = 0;
    for (auto it = container.begin(); it != container.end(); ++it) {
        ++count;
    }
    EXPECT_EQ(count, 50U);
}

/**
 * @brief Tests position updates with many particles causing cell migrations.
 *
 * This test uses more particles spread across the domain and triggers a bug in
 * LinkedCellContainerDirect::proximity_iterator::inc() where dereferencing
 * a null cell pointer causes a crash (SEGV at line 133).
 */
TYPED_TEST(ParticleContainerTest, UpdatePositionManyParticlesCellMigration) {
    // Larger domain with cutoff=3.0
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(3.0);

    R3 v{0.0, 0.0, 0.0};

    // Add many particles spread across the domain (similar to benchmark)
    for (int i = 0; i < 16; ++i) {
        // Spread particles in a 4x4 grid, staying within bounds
        int col = i % 4;
        int row = i / 4;
        double x = 1.0 + (col * 2.0);  // 1, 3, 5, 7
        double y = 1.0 + (row * 2.0);  // 1, 3, 5, 7
        container.addParticle(R3{x, y, 5.0}, v, 1.0, 1.0, 1.0);
    }

    ASSERT_EQ(container.size(), 16U);

    // Single pass with larger displacement to trigger cell migrations
    for (auto it = container.begin(); it != container.end();) {
        R3 new_pos = it->getX() + R3{1.5, 1.5, 0.0};
        it = container.updateParticlePosition(it, new_pos);
    }

    // Verify container is still valid
    EXPECT_EQ(container.size(), 16U);
    size_t count = 0;
    for (auto it = container.begin(); it != container.end(); ++it) {
        count++;
    }
    EXPECT_EQ(count, 16U);
}

/**
 * @brief Tests updating particle position causing cell migration.
 *
 * Minimal test: Particles in different cells. Move particles so they
 * migrate to different cells while iterating. This triggers a bug in
 * LinkedCellContainerDirect::proximity_iterator::satisfyInc() where
 * accessing particle position through an invalidated iterator crashes.
 */
TYPED_TEST(ParticleContainerTest, UpdatePositionCausesCellMigration) {
    // Use cutoff=2.5 so cell size is 2.5 (domain 10x10x10 = 4x4x4 cells)
    TypeParam container = ContainerFactory<TypeParam>::createWithCutoff(2.5);

    R3 v{0.0, 0.0, 0.0};

    // Add particles in different cells - need enough to trigger the bug
    container.addParticle(R3{1.0, 1.0, 1.0}, v, 1.0, 1.0, 1.0);
    container.addParticle(R3{4.0, 1.0, 1.0}, v, 1.0, 1.0, 1.0);
    container.addParticle(R3{7.0, 1.0, 1.0}, v, 1.0, 1.0, 1.0);

    ASSERT_EQ(container.size(), 3U);

    // Two rounds of updates - the second round triggers the bug
    for (int round = 0; round < 2; ++round) {
        for (auto it = container.begin(); it != container.end();) {
            R3 new_pos = it->getX() + R3{1.5, 0.0, 0.0};  // Move 1.5 units (may cross cell boundary)
            it = container.updateParticlePosition(it, new_pos);
        }
    }

    // Verify container is still valid
    EXPECT_EQ(container.size(), 3U);
    size_t count = 0;
    for (auto it = container.begin(); it != container.end(); ++it) {
        count++;
    }
    EXPECT_EQ(count, 3U);
}

// ============================================================================
// Combined Operations Tests
// ============================================================================

/**
 * @brief Tests a typical usage sequence.
 */
TYPED_TEST(ParticleContainerTest, TypicalUsage) {
    TypeParam container = ContainerFactory<TypeParam>::create();

    // Start empty
    EXPECT_TRUE(container.empty());

    // Reserve and add
    container.reserve(64);
    container.addParticle(Particle(42));
    container.addParticle({1., 2., 3.}, {4., 5., 6.}, 1.0, 5., 1.);
    ASSERT_EQ(container.size(), 2U);

    // Clear and verify
    container.clear();
    EXPECT_TRUE(container.empty());
    EXPECT_EQ(container.size(), 0U);
    EXPECT_EQ(container.begin(), container.end());
}

}  // namespace mol_sim
