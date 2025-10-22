#pragma once

#include <math.h>

#include <array>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

template <class T, size_t N>
  requires std::integral<T> || std::floating_point<T>
class Vector {
  std::array<T, N> _data;

 public:
  // constructors
  Vector() { _data.fill(T(0)); };
  Vector(std::array<T, N> init) : _data(init){};

  // cpy constr
  Vector(const Vector<T, N>& other) = default;
  // cpy assignment
  Vector<T, N>& operator=(const Vector<T, N>& other) = default;
  // mv constr
  Vector(Vector<T, N>&& other) = default;
  // mv assignment
  Vector<T, N>& operator=(Vector<T, N>&& other) = default;
  // dstr
  ~Vector() = default;

  // access
  T& operator[](size_t idx) { return _data[idx]; };
  const T& operator[](size_t idx) const { return _data[idx]; };
  T* data() { return _data.data(); }
  const T* data() const { return _data.data(); }

  // comparison
  bool operator==(const Vector<T, N>& other) const { return _data == other._data; };
  std::strong_ordering operator<=>(const Vector<T, N>& other) const { return euclidNorm() <=> other.euclidNorm(); };

  // point-wise arithmetic operators
  Vector<T, N> operator+(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] + other._data[i];
    return Vector<T, N>(new_array);
  };
  Vector<T, N> operator-(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] - other._data[i];
    return Vector<T, N>(new_array);
  };
  Vector<T, N> operator*(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] * other._data[i];
    return Vector<T, N>(new_array);
  };
  Vector<T, N> operator/(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] / other._data[i];
    return Vector<T, N>(new_array);
  };

  // other arithmetic operations
  Vector<T, N> operator*(const T scalar) const {  // Vector * Scalar
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = scalar * _data[i];
    return Vector<T, N>(new_array);
  }
  friend Vector<T, N> operator*(const T& s, const Vector<T, N>& v) { return v * s; }  // Scalar * Vector
  static T scalarProduct(const Vector& a, const Vector<T, N>& b) {
    return std::inner_product(a._data.begin(), a._data.end(), b._data.begin(), T(0));
  }
  T euclidNorm() const {
    return std::sqrt(std::accumulate(_data.begin(), _data.end(), T(0), [](auto a, auto b) { return a + b * b; }));
  }

  // conversions to arrays
  operator std::array<T, 3>() const { return _data; };
  std::array<T, 3>& operator*() { return _data; };
  const std::array<T, 3>& operator*() const { return _data; };

  // output
  std::string toString(const std::string& delimiter = ", ",
                       const std::array<std::string, 2>& surround = {"[", "]"}) const {
    std::stringstream strStream;
    strStream << surround[0];
    for (size_t i = 0; i < N; ++i) {
      if (i) strStream << ", ";
      strStream << _data[i];
    }
    strStream << surround[1];
    return strStream.str();
  };
};

template <class T>
  requires std::integral<T> || std::floating_point<T>
class Vector<T, 3> {
  std::array<T, 3> _data;

 public:
  // constructors
  Vector() : _data({T(0), T(0), T(0)}){};
  Vector(T x, T y, T z) : _data({x, y, z}){};
  Vector(T val) : _data({val, val, val}){};
  Vector(std::array<T, 3> init) : _data(init){};

  // cpy constr
  Vector(const Vector<T, 3>& other) = default;
  // cpy assignment
  Vector<T, 3>& operator=(const Vector<T, 3>& other) = default;
  // mv constr
  Vector(Vector<T, 3>&& other) = default;
  // mv assignment
  Vector<T, 3>& operator=(Vector<T, 3>&& other) = default;
  // dstr
  ~Vector() = default;

  // access
  T& operator[](size_t idx) { return _data[idx]; };
  const T& operator[](size_t idx) const { return _data[idx]; };
  T* data() { return _data.data(); }
  const T* data() const { return _data.data(); }

  T& x() { return _data[0]; };
  const T& x() const { return _data[0]; };
  T& y() { return _data[1]; };
  const T& y() const { return _data[1]; };
  T& z() { return _data[2]; };
  const T& z() const { return _data[2]; };

  // comparison
  bool operator==(const Vector<T, 3>& other) const { return _data == other._data; };
  std::strong_ordering operator<=>(const Vector<T, 3>& other) const { return euclidNorm() <=> other.euclidNorm(); };

  // point-wise arithmetic operators
  Vector<T, 3> operator+(const Vector<T, 3>& other) const {
    return Vector<T, 3>(_data[0] + other._data[0], _data[1] + other._data[1], _data[2] + other._data[2]);
  };
  Vector<T, 3> operator-(const Vector<T, 3>& other) const {
    return Vector<T, 3>(_data[0] - other._data[0], _data[1] - other._data[1], _data[2] - other._data[2]);
  };
  Vector<T, 3> operator*(const Vector<T, 3>& other) const {
    return Vector<T, 3>(_data[0] * other._data[0], _data[1] * other._data[1], _data[2] * other._data[2]);
  };
  Vector<T, 3> operator/(const Vector<T, 3>& other) const {
    return Vector<T, 3>(_data[0] / other._data[0], _data[1] / other._data[1], _data[2] / other._data[2]);
  };

  // other arithmetic operations
  Vector<T, 3> operator*(const T scalar) const {  // Vector * Scalar
    return Vector<T, 3>(_data[0] * scalar, _data[1] * scalar, _data[2] * scalar);
  };
  friend Vector<T, 3> operator*(const T& s, const Vector<T, 3>& v) { return v * s; }  // Scalar * Vector
  static T scalarProduct(const Vector<T, 3>& a, const Vector<T, 3>& b) {
    return (a._data[0] * b._data[0]) + (a._data[1] * b._data[1]) + (a._data[2] * b._data[2]);
  }
  T euclidNorm() const { return std::sqrt((_data[0] * _data[0]) + (_data[1] * _data[1]) + (_data[2] * _data[2])); };

  // conversion to array
  operator std::array<T, 3>() const { return _data; };
  std::array<T, 3>& operator*() { return _data; };
  const std::array<T, 3>& operator*() const { return _data; };

  // output
  std::string toString(const std::string& delimiter = ", ",
                       const std::array<std::string, 2>& surround = {"[", "]"}) const {
    std::stringstream strStream;
    strStream << surround[0] << x() << delimiter << y() << delimiter << z() << surround[1];
    return strStream.str();
  };
};

template <class T>
  requires std::integral<T> || std::floating_point<T>
class Vector<T, 2> {
  std::array<T, 2> _data;

 public:
  // constructors
  Vector() : _data({T(0), T(0)}){};
  Vector(T x, T y) : _data({x, y}){};
  Vector(T val) : _data({val, val}){};
  Vector(std::array<T, 2> init) : _data(init){};

  // cpy constr
  Vector(const Vector<T, 2>& other) = default;
  // cpy assignment
  Vector<T, 2>& operator=(const Vector<T, 2>& other) = default;
  // mv constr
  Vector(Vector<T, 2>&& other) = default;
  // mv assignment
  Vector<T, 2>& operator=(Vector<T, 2>&& other) = default;
  // dstr
  ~Vector() = default;

  // access
  T& operator[](size_t idx) { return _data[idx]; };
  const T& operator[](size_t idx) const { return _data[idx]; };

  T* data() { return _data.data(); }
  const T* data() const { return _data.data(); }

  T& x() { return _data[0]; };
  const T& x() const { return _data[0]; };
  T& y() { return _data[1]; };
  const T& y() const { return _data[1]; };

  // comparison
  bool operator==(const Vector<T, 2>& other) const { return _data == other._data; };
  std::strong_ordering operator<=>(const Vector<T, 2>& other) const { return euclidNorm() <=> other.euclidNorm(); };

  // point-wise arithmetic operators
  Vector<T, 2> operator+(const Vector<T, 2>& other) const {
    return Vector<T, 2>(_data[0] + other._data[0], _data[1] + other._data[1]);
  };
  Vector<T, 2> operator-(const Vector<T, 2>& other) const {
    return Vector<T, 2>(_data[0] - other._data[0], _data[1] - other._data[1]);
  };
  Vector<T, 2> operator*(const Vector<T, 2>& other) const {
    return Vector<T, 2>(_data[0] * other._data[0], _data[1] * other._data[1]);
  };
  Vector<T, 2> operator/(const Vector<T, 2>& other) const {
    return Vector<T, 2>(_data[0] / other._data[0], _data[1] / other._data[1]);
  };

  // other arithmetic operations
  Vector<T, 2> operator*(const T scalar) const {  // Vector * Scalar
    return Vector<T, 2>(_data[0] * scalar, _data[1] * scalar);
  };
  friend Vector<T, 2> operator*(const T& s, const Vector<T, 2>& v) { return v * s; }  // Scalar * Vector
  static T scalarProduct(const Vector<T, 2>& a, const Vector<T, 2>& b) {
    return a._data[0] * b._data[0] + a._data[1] * b._data[1];
  }
  T euclidNorm() const { return std::sqrt(_data[0] * _data[0] + _data[1] * _data[1]); };

  // conversion to arrays
  operator std::array<T, 2>() const { return _data; };
  std::array<T, 2>& operator*() { return _data; };
  const std::array<T, 2>& operator*() const { return _data; };

  // output
  std::string toString(const std::string& delimiter = ", ",
                       const std::array<std::string, 2>& surround = {"[", "]"}) const {
    std::stringstream strStream;
    strStream << surround[0] << x() << delimiter << y() << surround[1];
    return strStream.str();
  };
};

// output
template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& v) {
  os << v.toString();
  return os;
}

/**
 *
 * auto [x, y, z] = Vector<double, 3>()
 *
 * tuple_size , std::get
 *
 */
