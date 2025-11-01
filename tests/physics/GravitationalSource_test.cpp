#include <gtest/gtest.h>

#include "particles/Particle.h"
#include "physics/GravitationalForce.h"

namespace mol_sim {

class GravitationalForceTest : public testing::Test {
   protected:
    Particle p1;
    Particle p2;

    GravitationalForceTest() : p1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0), p2({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 2.0) {}
};

TEST_F(GravitationalForceTest, TestForceCalc) {
    R3 result = GravitationalForce().calculateForce(p1, p2);
    R3 expected = {2.0, 0.0, 0.0};
    EXPECT_EQ(result, expected);
}

}  // namespace mol_sim
