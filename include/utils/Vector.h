#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

namespace mol_sim {

/**
 * @brief Data type for Mathematical Vector
 *
 * Representation of a mathematical vector of dimension N and elements of type T (integers or floating point).
 * Enables access to all components (also via structured bindings), point-wise vector-vector operations, vector-scalar
 * (and scalar-vector) operations, other mathematical properties (euclidean norm, scalar product) and printing.
 *
 * @tparam T Component type, has to fulfill either  `std::integral` or `std::floating_point`
 * @tparam N Dimension of the Vector, meaning the number of components
 */
template <class T, size_t N>
    requires std::integral<T> || std::floating_point<T>
class Vector {
    /**
     * Data of the Vector
     */
    std::array<T, N> data_;  // NOLINT

   public:
    /**
     * @brief Default constructor produces Zero-Vector of required type
     */
    Vector() { data_.fill(T(0)); };

    /**
     * @brief Constructor copies given array into the components of a Vector
     *
     * @param init array of fitting type and size, copied to the data of the Vector
     */
    Vector(std::array<T, N> init) : data_(init) {};

    /**
     * @brief Constructor copies N given values into the components of a Vector
     *
     * @param inits arbitrary number (smaller or equal N) of elements or references of fitting type
     */
    template <typename... Ts>
        requires(... && (std::integral<std::decay_t<Ts>> || std::floating_point<std::decay_t<Ts>>))
    Vector(Ts&&... inits) : data_{inits...} {}

    // memory management

    /**
     * @brief Default Copy Constructor copies all components of a given Vector
     *
     * @param other Vector to copy
     */
    Vector(const Vector<T, N>& other) = default;

    /**
     * @brief Default Copy Assignment copies all components of a given Vector
     *
     * @param other Vector to copy
     */
    Vector<T, N>& operator=(const Vector<T, N>& other) = default;

    /**
     * @brief Default Move Constructor moves all components of a given Vector
     *
     * @param other Vector from which to move
     */
    Vector(Vector<T, N>&& other) = default;

    /**
     * @brief Default Move Assignment moves all components of a given Vector
     *
     * @param other Vector from which to move
     */
    Vector<T, N>& operator=(Vector<T, N>&& other) = default;

    /**
     * @brief Default Destructor destroys all components and underlying array
     */
    ~Vector() = default;

    // access

    T& operator[](size_t idx) { return data_[idx]; };
    const T& operator[](size_t idx) const { return data_[idx]; };
    /**
     * @brief Const pointer to the first element of the Vector, for usage as a C-style array
     *
     * @return Const pointer to the first element of the Vector
     */
    T* data() { return data_.data(); }
    /**
     * @brief Pointer to the first element of the Vector, for usage as a C-style array
     *
     * @return Pointer to the first element of the Vector
     */
    [[nodiscard]] const T* data() const { return data_.data(); }

    // comparison
    bool operator==(const Vector<T, N>& other) const { return data_ == other.data_; };
    std::strong_ordering operator<=>(const Vector<T, N>& other) const { return euclidNorm() <=> other.euclidNorm(); };

    // point-wise arithmetic operators

    /**
     * @brief Point-wise Addition of two Vectors
     *
     * @param other Vector of same type and dimension to add
     *
     * @return New Vector representing the sum of the two Vectors
     */
    Vector<T, N> operator+(const Vector<T, N>& other) const {
        Vector<T, N> result;
        for (size_t i = 0; i < N; i++) {
            result[i] = data_[i] + other.data_[i];
        }
        return result;
    };

    /**
     * @brief Point-wise Subtraction of two Vectors
     *
     * @param other Vector of same type and dimension to subtract
     *
     * @return New Vector representing the difference of the two Vectors
     */
    Vector<T, N> operator-(const Vector<T, N>& other) const {
        Vector<T, N> result;
        for (size_t i = 0; i < N; i++) {
            result[i] = data_[i] - other.data_[i];
        }
        return result;
    };

    /**
     * @brief Point-wise Multiplication of two Vectors (NOT INNER PRODUCT)
     *
     * @param other Vector of same type and dimension to multiply component-wise
     *
     * @return New Vector with the component-wise product of the two Vectors
     */
    Vector<T, N> operator*(const Vector<T, N>& other) const {
        Vector<T, N> result;
        for (size_t i = 0; i < N; i++) {
            result[i] = data_[i] * other.data_[i];
        }
        return result;
    };

    /**
     * @brief Point-wise Division of two Vectors
     *
     * @param other Vector of same type and dimension to devise by
     *
     * @return New Vector representing the quotient of the two Vectors
     */
    Vector<T, N> operator/(const Vector<T, N>& other) const {
        Vector<T, N> result;
        for (size_t i = 0; i < N; i++) {
            result[i] = data_[i] / other.data_[i];
        }
        return result;
    };

    /**
     * @brief Compound point-wise Addition of a Vector to this Vector
     *
     * @param other Vector of same type and dimension to add
     *
     * @return This Vector representing the sum of the two Vectors
     */
    Vector<T, N>& operator+=(const Vector<T, N>& other) {
        for (size_t i = 0; i < N; i++) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    /**
     * @brief Compound point-wise Subtraction of a Vector from this Vector
     *
     * @param other Vector of same type and dimension to subtract
     *
     * @return This Vector representing the difference of the two Vectors
     */
    Vector<T, N>& operator-=(const Vector<T, N>& other) {
        for (size_t i = 0; i < N; i++) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    // other arithmetic operations

    /**
     * @brief Scalar Multiplication (Vector * Scalar)
     *
     * @param scalar Scalar value of same type
     *
     * @return New Vector representing the scaled Vector
     */
    Vector<T, N> operator*(const T scalar) const {
        std::array<T, N> new_array;
        for (size_t i = 0; i < N; i++) {
            new_array[i] = scalar * data_[i];
        }
        return Vector<T, N>(new_array);
    }

    /**
     * @brief Scalar Multiplication (Scalar * Vector)
     *
     * @param s Scalar value of same type as Vector
     * @param v Vector
     *
     * @return New Vector representing the scaled Vector
     */
    friend Vector<T, N> operator*(const T& s, const Vector<T, N>& v) { return v * s; }

    /**
     * @brief Compound Scalar Multiplication (Vector *= Scalar)
     *
     * @param scalar Scalar value of same type
     *
     * @return This Vector representing the scaled Vector
     */
    Vector<T, N>& operator*=(const T scalar) {
        for (size_t i = 0; i < N; i++) {
            data_[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief Inner product of two Vectors < a | b >
     *
     * @param a Vector a on left-hand-side
     * @param b Vector b on right-hand-side
     *
     * @return Scalar number of type T
     */
    static T scalarProduct(const Vector& a, const Vector<T, N>& b) {
        return std::inner_product(a.data_.begin(), a.data_.end(), b.data_.begin(), T(0));
    }

    /**
     * @brief Euclidean Norm of the Vector (L2 norm)
     *
     * @return Scalar number of type T
     */
    [[nodiscard]] T euclidNorm() const {
        return std::sqrt(std::accumulate(data_.begin(), data_.end(), T(0), [](auto a, auto b) { return a + (b * b); }));
    }

    /**
     * @brief Squared Euclidean Norm of the Vector
     *
     * @return Scalar number of type T
     */
    [[nodiscard]] T sqrEuclidNorm() const {
        return std::accumulate(data_.begin(), data_.end(), T(0), [](auto a, auto b) { return a + (b * b); });
    }

    // conversions to arrays

    /**
     * @brief Implicit comversion to a std::array<T, 3>
     */
    operator std::array<T, 3>() const { return data_; };

    /**
     * @brief "Dereferencing" a Vector results in a std::array<T, 3>
     */
    std::array<T, 3>& operator*() { return data_; };

    /**
     * @brief "Dereferencing" a const Vector results in a const std::array<T, 3>
     */
    const std::array<T, 3>& operator*() const { return data_; };

    // output

    /**
     * @brief Conversion to a std::string with custom delimiter and brackets
     *
     * @param delimiter Custom delimiter between the components (default: `", "`)
     * @param surround 2-element array containing the surrounding strings of the Vector (default: `{"[", "]"}`)
     *
     * @return String representation of the Vector
     */
    [[nodiscard]] std::string toString(const std::string& delimiter = ", ",
                                       const std::array<std::string, 2>& surround = {"[", "]"}) const {
        std::stringstream str_stream;
        str_stream << surround[0];
        for (size_t i = 0; i < N; ++i) {
            if (i != 0U) {
                str_stream << delimiter;
            }
            str_stream << data_[i];
        }
        str_stream << surround[1];
        return str_stream.str();
    };
};

// output

/**
 * @brief Allows for printing with streams
 *
 * @param os Output stream
 * @param v Vector to print into the stream
 *
 * @return Output stream
 */
template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& v) {
    os << v.toString();
    return os;
}

}  // namespace mol_sim

// --------------------  STRUCTURED BINDINGS  --------------------

/**
 * @brief Allows for structured bindings
 * The following partial specializations of std templates allow for
 * auto [x, y, z] = vec;
 * with vec being a Vector<T, 3>
 */
namespace std {

template <class T, size_t N>
struct tuple_size<mol_sim::Vector<T, N>> : std::integral_constant<std::size_t, N> {
    static constexpr size_t VALUE = N;
};

template <size_t I, class T, size_t N>
struct tuple_element<I, mol_sim::Vector<T, N>> {
    using type = T;
};

}  // namespace std

namespace mol_sim {

template <size_t I, class T, size_t N>
constexpr T& get(Vector<T, N>& v) noexcept {
    return v[I];
}

template <size_t I, class T, size_t N>
constexpr const T& get(const Vector<T, N>& v) noexcept {
    return v[I];
}

template <size_t I, class T, size_t N>
constexpr T&& get(Vector<T, N>&& v) noexcept {
    return std::move(v[I]);
}

template <size_t I, class T, size_t N>
constexpr const T&& get(const Vector<T, N>&& v) noexcept {
    return std::move(v[I]);
}

/**
 * @brief Alias for a mathematical 3-dimensional real-valued (double-precision) vector.
 *
 * @see Vector
 */
using R3 = Vector<double, 3>;

}  // namespace mol_sim

// --------------------  STRUCTURED BINDINGS  --------------------

#endif
