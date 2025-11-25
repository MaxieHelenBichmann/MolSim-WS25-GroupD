#include "particles/container/LinkedCellContainer.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <set>
#include <unordered_set>
#include <vector>

#include "particles/ParticleContainer.h"
#include "utils/Vector.h"

namespace mol_sim {

static_assert(ParticleContainer<LinkedCellContainer>, "LinkedCellContainer must satisfy ParticleContainer concept");

/**
 * @brief Test Fixture for testing the LinkedCellContainer.
 *
 * Base Config/Data of the container is as follow:
 *
 */
class LinkedCellContainerTest : public testing::Test {
   protected:
    R3 domain_size;
    LinkedCellContainer particles_empty;
    LinkedCellContainer particles_full;
    Particle p0;
    Particle p1;
    Particle p2;
    Particle p3;
    double cutoff_radius_even = 2.5;

    LinkedCellContainerTest()
        : domain_size(R3{10.0, 10.0, 10.0}),
          particles_empty(LinkedCellContainer(domain_size, 2.5)),
          particles_full(LinkedCellContainer(domain_size, 2.5)),
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
TEST_F(LinkedCellContainerTest, testAccessSubscript) {
    EXPECT_TRUE(particles_full[0] == p0 || particles_full[0] == p1 || particles_full[0] == p2 ||
                particles_full[0] == p3);
    EXPECT_TRUE(particles_full[2] == p0 || particles_full[2] == p1 || particles_full[2] == p2 ||
                particles_full[2] == p3);
}

/**
 * @brief Tests correct behaviour of method size().
 */
TEST_F(LinkedCellContainerTest, testSize) {
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(0));
    particles_empty.addParticle(Particle(5));
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(4));
}

/**
 * @brief Tests correct behaviour of method empty().
 */
TEST_F(LinkedCellContainerTest, testEmpty) {
    EXPECT_FALSE(particles_full.empty());
    EXPECT_TRUE(particles_empty.empty());
    particles_empty.addParticle(Particle(4));
    EXPECT_FALSE(particles_empty.empty());
}

// modify

/**
 * @brief Tests correct behaviour of method clear().
 */
TEST_F(LinkedCellContainerTest, testClear) {
    ASSERT_EQ(particles_full.size(), static_cast<size_t>(4));
    particles_full.clear();
    EXPECT_TRUE(particles_full.empty());
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(0));
    EXPECT_EQ(particles_full.begin(), particles_full.end());
}

/**
 * @brief Tests correct behaviour of method reserve(size_t n).
 */
TEST_F(LinkedCellContainerTest, testReserve) {
    particles_empty.reserve(64);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(0));
}

/**
 * @brief Tests correct behaviour of the method overload addParticle(Particle&& value).
 */
TEST_F(LinkedCellContainerTest, testAddParticleRval) {
    Particle tmp(6);
    particles_empty.addParticle(std::move(tmp));
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method overload addParticle(const Particle& value).
 */
TEST_F(LinkedCellContainerTest, testAddParticleConstLval) {
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
TEST_F(LinkedCellContainerTest, testAddParticleEmplaceNoType) {
    particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 5., 1.);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method overload
 * addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type).
 */
TEST_F(LinkedCellContainerTest, testAddParticleEmplace) {
    particles_empty.addParticle({1., 2., 3.}, {4., 5., 6.}, 3.14, 8., 9., 69);
    EXPECT_EQ(particles_empty.size(), static_cast<size_t>(1));
}

/**
 * @brief Tests correct behaviour of the method eraseParticle(Particle* p).
 */
TEST_F(LinkedCellContainerTest, testEraseParticle) {
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(4));
    Particle* p_ptr = &particles_full[2];
    particles_full.eraseParticle(p_ptr);
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(3));

    Particle p(42);
    Particle* p_ptr_invalid = &p;
    particles_full.eraseParticle(p_ptr_invalid);
    EXPECT_EQ(particles_full.size(), static_cast<size_t>(3));
}

// normal iterators

/**
 * @brief Tests correct behaviour of required non-const iterator begin().
 */
TEST_F(LinkedCellContainerTest, testBeginIterator) {
    EXPECT_NE(particles_full.begin(), particles_full.end());
    EXPECT_EQ(particles_empty.begin(), particles_empty.end());
}

/**
 * @brief Tests correct behaviour of required const iterators begin() and cbegin().
 */
TEST_F(LinkedCellContainerTest, testBeginConstInterator) {
    // begin()
    const LinkedCellContainer& cc = particles_full;
    EXPECT_NE(cc.begin(), cc.end());

    // cbegin()
    EXPECT_NE(particles_full.cbegin(), particles_full.cend());
}

// proximity iterators

/**
 * @brief Tests correct behaviour of proximity iterator with infinite radius.
 */
TEST_F(LinkedCellContainerTest, testProximityIteratorInfiniteRadius) {
    R3 v{0.0, 0.0, 0.0};
    particles_empty.addParticle(R3{2.0, 2.0, 2.0}, v, 1.0, 1.0, 1.0);
    particles_empty.addParticle(R3{4.0, 4.0, 4.0}, v, 1.0, 1.0, 1.0);
    particles_empty.addParticle(R3{6.0, 6.0, 6.0}, v, 1.0, 1.0, 1.0);

    R3 center{3.0, 3.0, 3.0};
    double radius = std::numeric_limits<double>::infinity();

    auto it = particles_empty.proximityBegin(center, radius);
    auto end = particles_empty.proximityEnd(center, radius);

    size_t count = 0;
    while (it != end) {
        EXPECT_LE((it->getX() - center).euclidNorm(), radius);
        ++it;
        ++count;
    }
    EXPECT_EQ(count, 3);  // Assuming only three particles are within the radius
}

/**
 * @brief Tests correct behaviour of proximity iterator with finite radius.
 */
TEST_F(LinkedCellContainerTest, testProximityIterator) {
    R3 v{0.0, 0.0, 0.0};
    particles_empty.addParticle(R3{2.5, 3.0, 3.0}, v, 1.0, 1.0, 1.0);
    particles_empty.addParticle(R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0);
    particles_empty.addParticle(R3{6.0, 6.0, 6.0}, v, 1.0, 1.0, 1.0);

    R3 center{3.0, 3.0, 3.0};
    double radius = 1.0;

    auto it = particles_empty.proximityBegin(center, radius);
    auto end = particles_empty.proximityEnd(center, radius);

    size_t count = 0;
    while (it != end) {
        EXPECT_LE((it->getX() - center).euclidNorm(), radius);
        ++it;
        ++count;
    }
    EXPECT_EQ(count, 1);  // Assuming only one particle is within the radius
}

// ParticleContainer: complex tests

/**
 * @brief Tests a sequence of read and write operations on a LinkedCellContainer.
 */
TEST(LinkedCellContainer, testSizeEmptyClearReserve) {
    LinkedCellContainer c({15., 10., 10.}, 2.);
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

// LinkedCellContainer: specific tests

/**
 * @brief Tests correct behaviour of method fitsDomain.
 */
TEST_F(LinkedCellContainerTest, testFitDomain) {
    R3 inside{5.0, 5.0, 5.0};
    EXPECT_TRUE(particles_empty.fitsDomain(inside));

    R3 far_outside{1000.0, 0.0, 0.0};
    EXPECT_FALSE(particles_empty.fitsDomain(far_outside));

    R3 close_inside{9.7, 10.0, 10.0};
    EXPECT_TRUE(particles_empty.fitsDomain(close_inside));

    R3 close_outside{10.0, 10.5, 10.0};
    EXPECT_FALSE(particles_empty.fitsDomain(close_outside));

    R3 on_min{0.0, 0.0, 0.0};
    R3 on_max{10.0, 10.0, 10.0};
    EXPECT_TRUE(particles_empty.fitsDomain(on_min));
    EXPECT_TRUE(particles_empty.fitsDomain(on_max));

    R3 outside_x{-0.1, 5.0, 5.0};
    R3 outside_y{5.0, -0.1, 5.0};
    R3 outside_z{5.0, 5.0, 10.1};
    EXPECT_FALSE(particles_empty.fitsDomain(outside_x));
    EXPECT_FALSE(particles_empty.fitsDomain(outside_y));
    EXPECT_FALSE(particles_empty.fitsDomain(outside_z));
}

/**
 * @brief Tests correct behaviour of method fitsContainer, which also considers halo cells.
 */
TEST_F(LinkedCellContainerTest, testFitContainer) {
    R3 inside{5.0, 5.0, 5.0};
    EXPECT_TRUE(particles_empty.fitsContainer(inside));

    R3 far_outside{1000.0, 0.0, 0.0};
    EXPECT_FALSE(particles_empty.fitsContainer(far_outside));

    R3 halo{-0.5, 5.0, 5.0};
    EXPECT_TRUE(particles_empty.fitsContainer(halo));

    R3 close_inside{-2.5, -2.5, -2.5};
    EXPECT_TRUE(particles_empty.fitsContainer(close_inside));

    R3 close_outside{12.5, 12.5, 12.7};
    EXPECT_FALSE(particles_empty.fitsContainer(close_outside));
}

/**
 * @brief Tests correct behaviour of the method updateParticlePosition, with different scenarios.
 */
TEST_F(LinkedCellContainerTest, testUpdateParticlePosition) {
    R3 v{0.0, 0.0, 0.0};

    Particle p_inside{R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle p_move_outside{R3{9.0, 9.0, 9.0}, v, 1.0, 1.0, 1.0};
    Particle p_move_halo{R3{1.0, 1.0, 1.0}, v, 1.0, 1.0, 1.0};

    particles_empty.addParticle(p_inside);
    particles_empty.addParticle(p_move_outside);
    particles_empty.addParticle(p_move_halo);

    ASSERT_EQ(particles_empty.size(), 3);

    // Move inside -> inside
    particles_empty.updateParticlePosition(particles_empty.begin(), R3{6.0, 6.0, 6.0});
    EXPECT_TRUE(particles_empty.begin()->getX() == (R3{6.0, 6.0, 6.0}));

    // Move inside -> outside
    particles_empty.updateParticlePosition(particles_empty.begin() + 1, R3{15.0, 15.0, 15.0});
    EXPECT_EQ(particles_empty.size(), 2);

    // Move inside -> halo
    particles_empty.updateParticlePosition(particles_empty.begin() + 1, R3{-1.0, -1.0, -1.0});
    EXPECT_TRUE((particles_empty.begin() + 1)->getX() == (R3{-1.0, -1.0, -1.0}));
}

/**
 * @brief Tests correct behaviour of the boundary iterator for all and a specific side.
 */
TEST_F(LinkedCellContainerTest, testBoundaryIterator) {  // NOLINT
    R3 v{0.0, 0.0, 0.0};

    Particle par_inside{R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_boundary_left{R3{0.1, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_boundary_right{R3{9.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};

    particles_empty.addParticle(par_inside);
    particles_empty.addParticle(par_boundary_left);
    particles_empty.addParticle(par_boundary_right);

    ASSERT_EQ(particles_empty.size(), 3);

    // All boundary sides
    std::vector<Particle> boundary_all;
    auto it = particles_empty.boundaryBegin();
    auto end = particles_empty.boundaryEnd();
    for (; it != end; ++it) {
        boundary_all.push_back(*it);
    }

    EXPECT_EQ(boundary_all.size(), 2);
    EXPECT_TRUE(std::count(boundary_all.begin(), boundary_all.end(), par_boundary_left) > 0);
    EXPECT_TRUE(std::count(boundary_all.begin(), boundary_all.end(), par_boundary_right) > 0);
    EXPECT_FALSE(std::count(boundary_all.begin(), boundary_all.end(), par_inside) > 0);

    // Only RIGHT boundary
    std::vector<Particle> boundary_right;
    for (auto it = particles_empty.boundaryBegin({BoundaryLocation::RIGHT});
         it != particles_empty.boundaryEnd({BoundaryLocation::RIGHT}); ++it) {
        boundary_right.push_back(*it);
    }

    EXPECT_EQ(boundary_right.size(), 1);
    EXPECT_FALSE(std::count(boundary_right.begin(), boundary_right.end(), par_boundary_left) > 0);
    EXPECT_TRUE(std::count(boundary_right.begin(), boundary_right.end(), par_boundary_right) > 0);
    EXPECT_FALSE(std::count(boundary_right.begin(), boundary_right.end(), par_inside) > 0);
}

/**
 * @brief Tests correct behaviour of the halo iterator for all and a specific side.
 */
TEST_F(LinkedCellContainerTest, testHaloIterator) {  // NOLINT
    R3 v{0.0, 0.0, 0.0};

    Particle par_inside{R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_boundary{R3{2.0, 2.0, 2.0}, v, 1.0, 1.0, 1.0};
    Particle par_halo_left{R3{-0.1, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_halo_right{R3{12.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};

    particles_empty.addParticle(par_inside);
    particles_empty.addParticle(par_boundary);
    particles_empty.addParticle(par_halo_left);
    particles_empty.addParticle(par_halo_right);

    ASSERT_EQ(particles_empty.size(), 4);

    // All boundary sides
    std::vector<Particle> halo_all;
    auto it = particles_empty.haloBegin();
    auto end = particles_empty.haloEnd();
    for (; it != end; ++it) {
        halo_all.push_back(*it);
    }

    EXPECT_EQ(halo_all.size(), 2);
    EXPECT_TRUE(std::count(halo_all.begin(), halo_all.end(), par_halo_left) > 0);
    EXPECT_TRUE(std::count(halo_all.begin(), halo_all.end(), par_halo_right) > 0);
    EXPECT_FALSE(std::count(halo_all.begin(), halo_all.end(), par_inside) > 0);
    EXPECT_FALSE(std::count(halo_all.begin(), halo_all.end(), par_boundary) > 0);
    // Only RIGHT boundary
    std::vector<Particle> halo_right;
    for (auto it = particles_empty.haloBegin({BoundaryLocation::RIGHT});
         it != particles_empty.haloEnd({BoundaryLocation::RIGHT}); ++it) {
        halo_right.push_back(*it);
    }

    EXPECT_EQ(halo_right.size(), 1);
    EXPECT_FALSE(std::count(halo_right.begin(), halo_right.end(), par_halo_left) > 0);
    EXPECT_TRUE(std::count(halo_right.begin(), halo_right.end(), par_halo_right) > 0);
    EXPECT_FALSE(std::count(halo_right.begin(), halo_right.end(), par_inside) > 0);
    EXPECT_FALSE(std::count(halo_right.begin(), halo_right.end(), par_boundary) > 0);
}

}  // namespace mol_sim
