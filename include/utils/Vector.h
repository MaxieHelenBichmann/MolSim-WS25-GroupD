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
 */
template <class T, size_t N>
  requires std::integral<T> || std::floating_point<T>
class Vector {
  /**
   * Data of the Vector
   */
  
  std::array<T, N> _data; //NOLINT

 public:
  /**
   * @brief Default constructor produces Zero-Vector of required type
   */
  Vector() { _data.fill(T(0)); };

  /**
   * @brief Constructor copies given array into the components of a Vector
   */
  Vector(std::array<T, N> init) : _data(init) {};

  /**
   * @brief Constructor copies N given values into the components of a Vector
   */
  template <typename... Ts>
    requires(... && (std::integral<Ts> || std::floating_point<Ts>))
  Vector(Ts&&... inits) : _data{inits...} {}

  // memory management

  /**
   * @brief Default Copy Constructor copies all components of a given Vector
   */
  Vector(const Vector<T, N>& other) = default;

  /**
   * @brief Default Copy Assignment copies all components of a given Vector
   */
  Vector<T, N>& operator=(const Vector<T, N>& other) = default;

  /**
   * @brief Default Move Constructor moves all components of a given Vector
   */
  Vector(Vector<T, N>&& other) = default;

  /**
   * @brief Default Move Assignment moves all components of a given Vector
   */
  Vector<T, N>& operator=(Vector<T, N>&& other) = default;

  /**
   * @brief Default Destructor destroys all components and underlying array
   */
  ~Vector() = default;

  // access

  T& operator[](size_t idx) { return _data[idx]; };
  const T& operator[](size_t idx) const { return _data[idx]; };
  T* data() { return _data.data(); }
  [[nodiscard]] const T* data() const { return _data.data(); }

  // comparison
  bool operator==(const Vector<T, N>& other) const { return _data == other._data; };
  std::strong_ordering operator<=>(const Vector<T, N>& other) const { return euclidNorm() <=> other.euclidNorm(); };

  // point-wise arithmetic operators

  /**
   * @brief Point-wise Addition of two Vectors
   */
  Vector<T, N> operator+(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] + other._data[i];
    }
    return Vector<T, N>(new_array);
  };

  /**
   * @brief Point-wise Subtraction of two Vectors
   */
  Vector<T, N> operator-(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] - other._data[i];
    }
    return Vector<T, N>(new_array);
  };

  /**
   * @brief Point-wise Multiplication of two Vectors (NOT INNER PRODUCT)
   */
  Vector<T, N> operator*(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] * other._data[i];
    }
    return Vector<T, N>(new_array);
  };

  /**
   * @brief Point-wise Division of two Vectors
   */
  Vector<T, N> operator/(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] / other._data[i];
    }
    return Vector<T, N>(new_array);
  };

  // other arithmetic operations

  /**
   * @brief Scalar Multiplication (Vector * Scalar)
   */
  Vector<T, N> operator*(const T scalar) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = scalar * _data[i];
    }
    return Vector<T, N>(new_array);
  }

  /**
   * @brief Scalar Multiplication (Scalar * Vector)
   */
  friend Vector<T, N> operator*(const T& s, const Vector<T, N>& v) { return v * s; }

  /**
   * @brief Inner product of two Vectors < a | b >
   */
  static T scalarProduct(const Vector& a, const Vector<T, N>& b) {
    return std::inner_product(a._data.begin(), a._data.end(), b._data.begin(), T(0));
  }

  /**
   * @brief Euclidean Norm of a Vector (L2 norm)
   */
  [[nodiscard]] T euclidNorm() const {
    return std::sqrt(std::accumulate(_data.begin(), _data.end(), T(0), [](auto a, auto b) { return a + (b * b); }));
  }

  // conversions to arrays

  /**
   * @brief Implicit comversion to a std::array<T, 3>
   */
  operator std::array<T, 3>() const { return _data; };

  /**
   * @brief "Dereferencing" a Vector results in a std::array<T, 3>
   */
  std::array<T, 3>& operator*() { return _data; };

  /**
   * @brief "Dereferencing" a const Vector results in a const std::array<T, 3>
   */
  const std::array<T, 3>& operator*() const { return _data; };

  // output

  /**
   * @brief Conversion to a std::string with custom delimiter and brackets
   */
  [[nodiscard]] std::string toString(const std::string& delimiter = ", ",
                                     const std::array<std::string, 2>& surround = {"[", "]"}) const {
    std::stringstream str_stream;
    str_stream << surround[0];
    for (size_t i = 0; i < N; ++i) {
      if (i != 0U) {
        str_stream << ", ";
      }
      str_stream << _data[i];
    }
    str_stream << surround[1];
    return str_stream.str();
  };
};

// output

/**
 * @brief Allows for printing with streams
 */
template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& v) {
  os << v.toString();
  return os;
}

}  // namespace mol_sim

/**
 * @brief Allows for structured bindings
 * The following partial specializations of std templates allow for
 * auto [x, y, z] = vec;
 * with vec being a Vector<T, 3>
 */
namespace std {

template <class T, size_t N>
struct tuple_size<mol_sim::Vector<T, N>> {
  static constexpr size_t VALUE = N;
};

template <size_t I, class T, size_t N>
constexpr T& get(mol_sim::Vector<T, N>&& v) {
  return v[I];
}

template <size_t I, class T, size_t N>
constexpr T get(const mol_sim::Vector<T, N>& v) {
  return v[I];
}

template <size_t I, class T, size_t N>
struct tuple_element<I, mol_sim::Vector<T, N>> {
  using type = T;
};

}  // namespace std

#endif
