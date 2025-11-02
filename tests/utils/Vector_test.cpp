#include "utils/Vector.h"

#include <gtest/gtest.h>

#include <compare>
#include <limits>

namespace mol_sim {

// Copy & Move semantics

/**
 * @brief Tests different constructors of Vector.
 */
TEST(VectorTest, TestDefaultConstructor) {
    // default constructor
    Vector<double, 10> r10 = Vector<double, 10>();
    Vector<int, 4> z4 = Vector<int, 4>();

    for (size_t idx = 0; idx < 10; idx++) {
        EXPECT_DOUBLE_EQ(r10[idx], 0.);
    }
    for (size_t idx = 0; idx < 4; idx++) {
        EXPECT_EQ(z4[idx], 0);
    }
}

TEST(VectorTest, TestArrayConstructor) {
    // array init
    std::array<double, 3> arr{1.5, -2.0, 3.25};
    Vector<double, 3> r3(arr);

    for (size_t idx = 0; idx < 3; idx++) {
        EXPECT_DOUBLE_EQ(r3[idx], arr[idx]);
    }
}

TEST(VectorTest, TestConstructorArbitraryValues) {
    // arbitrary argument list
    Vector<int, 10> z10 = Vector<int, 10>(20502, 1301, 5, 0, 24, 66, 69, 74, 9, 10);
    EXPECT_EQ(z10[0], 20502);
    EXPECT_EQ(z10[1], 1301);
    EXPECT_EQ(z10[2], 5);
    EXPECT_EQ(z10[9], 10);

    Vector<int, 5> z5 = Vector<int, 5>(69, 1337);
    EXPECT_EQ(z5[0], 69);
    EXPECT_EQ(z5[1], 1337);
    EXPECT_EQ(z5[2], 0);
    EXPECT_EQ(z5[3], 0);
    EXPECT_EQ(z5[4], 0);

    Vector<double, 4> r4 = Vector<double, 4>(3., 4.);
    EXPECT_DOUBLE_EQ(r4[0], 3.);
    EXPECT_DOUBLE_EQ(r4[1], 4.);
    EXPECT_DOUBLE_EQ(r4[2], 0.);
    EXPECT_DOUBLE_EQ(r4[3], 0.);
}

/**
 * @brief Tests copy constructor and assignment of Vector.
 */
TEST(VectorTest, TestCopy) {
    // copy constructor
    Vector<int, 3> a(1, 2, 3);
    Vector<int, 3> b = a;
    EXPECT_TRUE(a == b);

    // copy assignment
    Vector<int, 3> c;
    c = a;
    EXPECT_TRUE(a == c);
}

/**
 * @brief Tests move constructor and assignment of Vector.
 */
TEST(VectorTest, TestMove) {
    Vector<double, 3> r3(4.7, 5., 6.9);
    Vector<int, 3> z3(4, 5, 6);

    // move constructor
    Vector<double, 3> moved(r3);
    EXPECT_DOUBLE_EQ(moved[0], 4.7);
    EXPECT_DOUBLE_EQ(moved[1], 5.);
    EXPECT_DOUBLE_EQ(moved[2], 6.9);

    // move assignment
    Vector<int, 3> dst;
    dst = z3;
    EXPECT_EQ(dst[0], 4);
    EXPECT_EQ(dst[1], 5);
    EXPECT_EQ(dst[2], 6);
}

// Access

/**
 * @brief Tests correct access with the subscript operator.
 */
TEST(VectorTest, TestAccessSubscript) {
    Vector<double, 10> r10(1.125, 0.66, 12.3, -45., 1337., 42.8, -8.1, NAN, 66.6, 0.22);
    Vector<double, 3> r3(13.5, 69., .01);
    Vector<int, 3> z3(700, 44, 222);

    EXPECT_DOUBLE_EQ(r3[0], 13.5);
    EXPECT_DOUBLE_EQ(r3[1], 69.);
    EXPECT_DOUBLE_EQ(r3[2], .01);

    EXPECT_DOUBLE_EQ(r10[1], 0.66);
    EXPECT_DOUBLE_EQ(r10[6], -8.1);

    EXPECT_EQ(z3[0], 700);
    EXPECT_EQ(z3[1], 44);
    EXPECT_EQ(z3[2], 222);
}

/**
 * @brief Tests correct access to components using structured bindings.
 *
 */
TEST(VectorTest, TestAccessStructuresBindings) {
    Vector<double, 3> r3(6.1, -2., 50.);
    Vector<int, 2> r2(4, 3);

    auto [x, y, z] = r3;
    EXPECT_DOUBLE_EQ(x, 6.1);
    EXPECT_DOUBLE_EQ(y, -2.);
    EXPECT_DOUBLE_EQ(z, 50.);

    auto [a, b] = r2;
    EXPECT_EQ(a, 4);
    EXPECT_EQ(b, 3);
}

// Comparison

/**
 * @brief Tests equality and ordering of different Vectors.
 */
TEST(VectorTest, TestEquality) {
    Vector<int, 3> a(1, 2, 3);
    Vector<int, 3> b(1, 2, 3);
    Vector<int, 3> c(1, 2, 4);

    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

/**
 * @brief Tests ordering of different Vectors.
 */
TEST(VectorTest, TestOrdering) {
    Vector<int, 2> a(3, 4);
    Vector<int, 2> b(6, 8);
    Vector<int, 2> c(4, 3);

    EXPECT_EQ((a <=> b), std::strong_ordering::less);
    EXPECT_EQ((b <=> a), std::strong_ordering::greater);
    EXPECT_EQ((a <=> c), std::strong_ordering::equal);
}

// Arithmetic Operations

/**
 * @brief Tests correctness of point-wise addition of two Vectors.
 */
TEST(VectorTest, TestPointwiseAddition) {
    Vector<double, 3> a_r3(13.5, 69., .01);
    Vector<double, 3> b_r3(6.1, -2., 50.);

    Vector<int, 2> a_z2(22, 61);
    Vector<int, 2> b_z2(74, 9);

    EXPECT_DOUBLE_EQ((a_r3 + b_r3)[0], 19.6);
    EXPECT_DOUBLE_EQ((a_r3 + b_r3)[1], 67.);
    EXPECT_DOUBLE_EQ((a_r3 + b_r3)[2], 50.01);

    EXPECT_EQ((a_z2 + b_z2)[0], 96);
    EXPECT_EQ((a_z2 + b_z2)[1], 70);

    EXPECT_TRUE((a_z2 + b_z2) == (b_z2 + a_z2));
}

/**
 * @brief Tests correctness of point-wise subtraction of two Vectors.
 */
TEST(VectorTest, TestPointwiseSubtraction) {
    Vector<double, 3> a_r3(13.5, 69., .01);
    Vector<double, 3> b_r3(6.1, -2., 50.);

    Vector<int, 2> a_z2(22, 61);
    Vector<int, 2> b_z2(74, 9);

    EXPECT_DOUBLE_EQ((a_r3 - b_r3)[0], 7.4);
    EXPECT_DOUBLE_EQ((a_r3 - b_r3)[1], 71.);
    EXPECT_DOUBLE_EQ((a_r3 - b_r3)[2], -49.99);

    EXPECT_EQ((a_z2 - b_z2)[0], -52);
    EXPECT_EQ((a_z2 - b_z2)[1], 52);
}

/**
 * @brief Tests correctness of point-wise multiplication of two Vectors.
 */
TEST(VectorTest, TestPointwiseMultiplication) {
    Vector<double, 3> a_r3(13.5, 69., .01);
    Vector<double, 3> b_r3(6.1, -2., 50.);

    Vector<int, 2> a_z2(22, 61);
    Vector<int, 2> b_z2(74, 9);

    EXPECT_DOUBLE_EQ((a_r3 * b_r3)[0], 82.35);
    EXPECT_DOUBLE_EQ((a_r3 * b_r3)[1], -138.);
    EXPECT_DOUBLE_EQ((a_r3 * b_r3)[2], 0.5);

    EXPECT_EQ((a_z2 * b_z2)[0], 1628);
    EXPECT_EQ((a_z2 * b_z2)[1], 549);

    EXPECT_TRUE((a_z2 * b_z2) == (b_z2 * a_z2));
}

/**
 * @brief Tests correctness of point-wise division of two Vectors.
 */
TEST(VectorTest, TestPointwiseDivision) {
    Vector<double, 3> a_r3(13.5, 69., .01);
    Vector<double, 3> b_r3(0.5, -2., 50.);

    Vector<int, 2> a_z2(22, 61);
    Vector<int, 2> b_z2(74, 9);

    EXPECT_DOUBLE_EQ((a_r3 / b_r3)[0], 27.);
    EXPECT_DOUBLE_EQ((a_r3 / b_r3)[1], -34.5);
    EXPECT_DOUBLE_EQ((a_r3 / b_r3)[2], 0.0002);

    EXPECT_EQ((a_z2 / b_z2)[0], 0);
    EXPECT_EQ((a_z2 / b_z2)[1], 6);
}

/**
 * @brief Tests correctness of scalar multiplication of two Vectors.
 */
TEST(VectorTest, TestInnerProduct) {
    Vector<int, 3> a_z3(1, 2, 3);
    Vector<int, 3> b_z3(4, -1, 0);
    Vector<double, 2> a_r2(4.5, 1.125);
    Vector<double, 2> b_r2(0.01, 56.7);

    EXPECT_EQ((Vector<int, 3>::scalarProduct(a_z3, b_z3)), 2);
    EXPECT_EQ((Vector<double, 2>::scalarProduct(a_r2, b_r2)), 63.8325);
}

/**
 * @brief Tests correctness of the scalar product of a Vector and a scalar in both orders.
 */
TEST(VectorTest, TestScalarMultiplication) {
    Vector<double, 3> r3(1.0, -2.0, 0.5);
    Vector<double, 10> r10(1.125, 0.66, 12.3, -45., 1337., 42.8, -8.1, 0., 66.6, 0.22);

    auto a = r3 * 2.0;
    EXPECT_DOUBLE_EQ(a[0], 2.0);
    EXPECT_DOUBLE_EQ(a[1], -4.0);
    EXPECT_DOUBLE_EQ(a[2], 1.0);

    auto b = 3.0 * r3;
    EXPECT_DOUBLE_EQ(b[0], 3.0);
    EXPECT_DOUBLE_EQ(b[1], -6.0);
    EXPECT_DOUBLE_EQ(b[2], 1.5);

    EXPECT_TRUE((5.5 * r3) == (r3 * 5.5));
}

/**
 * @brief Tests correctness of the euclidean norm of a Vector.
 */
TEST(VectorTest, TestEuclideanNorm) {
    Vector<double, 3> r3(3.0, 4.0, 12.0);
    Vector<int, 2> z2(3, 4);

    EXPECT_DOUBLE_EQ(r3.euclidNorm(), 13.0);
    EXPECT_EQ(z2.euclidNorm(), 5);
}

// Conversions

/**
 * @brief Tests implicit conversion of a Vector to a std::array.
 */
TEST(VectorTest, TestImplicitConversion) {
    Vector<int, 3> z3(7, 8, 9);
    std::array<int, 3> arr = z3;
    EXPECT_EQ(arr[0], 7);
    EXPECT_EQ(arr[1], 8);
    EXPECT_EQ(arr[2], 9);
}

/**
 * @brief Tests explicit conversion of a Vector to a C-style array.
 */
TEST(VectorTest, TestExplicitDeref) {
    // non-const
    Vector<int, 3> z3(1, 2, 3);
    std::array<int, 3>& arr = *z3;
    EXPECT_EQ(z3[1], 2);
    arr[1] = 42;
    EXPECT_EQ(z3[1], 42);

    // const
    const Vector<int, 3> c_z3(69, 1213);
    const std::array<int, 3>& carr = *c_z3;
    EXPECT_EQ(carr[0], 69);
    EXPECT_EQ(carr[1], 1213);
    EXPECT_EQ(carr[2], 0);
}

// Strings & Stream output

/**
 * @brief Tests correctness of the default string representation of a Vector.
 */
TEST(VectorTest, TestStringRepresentation) {
    Vector<int, 3> z3(1, 2, 3);
    Vector<double, 5> r5(3.0, 4.5, NAN, 5.e-2, std::numeric_limits<double>::infinity());

    EXPECT_STREQ(z3.toString().c_str(), "[1, 2, 3]");
    EXPECT_STREQ(r5.toString().c_str(), "[3, 4.5, nan, 0.05, inf]");
}

/**
 * @brief Tests correctness of a custom string representation of a Vector.
 */
TEST(VectorTest, TestCustomStringRepresentation) {
    Vector<int, 3> z3(1, 2, 3);
    Vector<double, 5> r5(3.9, 5.7);

    EXPECT_STREQ(z3.toString(";", {"{", "}"}).c_str(), "{1; 2; 3}");
    EXPECT_STREQ(r5.toString(" |").c_str(), "[3.9 | 5.7 | 0 | 0 | 0]");
}

/**
 * @brief Tests behavior of Vector when inserted in OStream.
 */
TEST(VectorTest, TestOStreamInsertion) {
    Vector<double, 3> r3(1.25, -0.5, 2.0);
    std::ostringstream oss;
    oss << r3;
    EXPECT_EQ(oss.str(), "[1.25, -0.5, 2]");
}

}  // namespace mol_sim
