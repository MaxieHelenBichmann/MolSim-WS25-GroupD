#pragma once

#include <math.h>

#include <array>
#include <iostream>
#include <numeric>

template <class T, size_t N>
  requires std::integral<T> || std::floating_point<T>
class Vector {
  std::array<T, N> _data;

  template <class U, std::size_t M>
  friend std::ostream& operator<<(std::ostream& os, const Vector<U, M>& v);

 public:
  // constructors
  Vector() { _data.fill(T(0)); };
  Vector(std::array<T, N> init) : _data(init) {};

  // access
  T& operator[](size_t idx) { return _data[idx]; };
  const T& operator[](size_t idx) const { return _data[idx]; };

  T* data() { return _data.data(); }
  const T* data() const { return _data.data(); }

  // point-wise arithmetic operators
  Vector operator+(const Vector& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] + other._data[i];
    return Vector<T, N>(new_array);
  };
  Vector operator-(const Vector& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] - other._data[i];
    return Vector<T, N>(new_array);
  };
  Vector operator*(const Vector& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] * other._data[i];
    return Vector<T, N>(new_array);
  };
  Vector operator/(const Vector& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) new_array[i] = _data[i] / other._data[i];
    return Vector<T, N>(new_array);
  };

  // other arithmetic operations
  static T scalarProduct(const Vector& a, const Vector& b) {
    return std::inner_product(a._data.begin(), a._data.end(), b._data.begin(), T(0));
  }
  T euclidNorm() const {
    return std::sqrt(std::accumulate(_data.begin(), _data.end(), T(0), [](auto a, auto b) { return a + b * b; }));
  }

  bool operator==(const Vector& other) const { return _data == other._data; };
  std::strong_ordering operator<=>(const Vector& other) const { return euclidNorm() <=> other.euclidNorm(); };

  operator std::array<T, 3>() const { return _data; };
  std::array<T, 3>& operator*() { return _data; };
  const std::array<T, 3>& operator*() const { return _data; };

  std::ostream& operator<<(std::ostream& os) const {
    os << "[";
    for (size_t i = 0; i < N; ++i) {
      if (i) os << ", ";
      os << _data[i];
    }
    return os << "]";
  };
};

template <class T>
  requires std::integral<T> || std::floating_point<T>
class Vector<T, 3> {
  std::array<T, 3> _data;

 public:
  // constructors
  Vector() : _data({T(0), T(0), T(0)}) {};
  Vector(T x, T y, T z) : _data({x, y, z}) {};
  Vector(T val) : _data({val, val, val}) {};
  Vector(std::array<T, 3> init) : _data(init) {};

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
  static T scalarProduct(const Vector<T, 3>& a, const Vector<T, 3>& b) {
    return a._data[0] * b._data[0] + a._data[1] * b._data[1] + a._data[2] * b._data[2];
  }
  T euclidNorm() const { return std::sqrt(_data[0] * _data[0] + _data[1] * _data[1] + _data[2] * _data[2]); };

  bool operator==(const Vector<T, 3>& other) const { return _data == other._data; };
  std::strong_ordering operator<=>(const Vector<T, 3>& other) const { return euclidNorm() <=> other.euclidNorm(); };

  operator std::array<T, 3>() const { return _data; };
  std::array<T, 3>& operator*() { return _data; };
  const std::array<T, 3>& operator*() const { return _data; };

  std::ostream& operator<<(std::ostream& os) const {
    os << "[" << _data[0] << ", " << _data[1] << ", " << _data[2] << "]";
    return os;
  }
};

template <class T>
  requires std::integral<T> || std::floating_point<T>
class Vector<T, 2> {
  std::array<T, 2> _data;

 public:
  // constructors
  Vector() : _data({T(0), T(0)}) {};
  Vector(T x, T y) : _data({x, y}) {};
  Vector(T val) : _data({val, val}) {};
  Vector(std::array<T, 2> init) : _data(init) {};

  // access
  T& operator[](size_t idx) { return _data[idx]; };
  const T& operator[](size_t idx) const { return _data[idx]; };

  T* data() { return _data.data(); }
  const T* data() const { return _data.data(); }

  T& x() { return _data[0]; };
  const T& x() const { return _data[0]; };
  T& y() { return _data[1]; };
  const T& y() const { return _data[1]; };

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
  static T scalarProduct(const Vector<T, 2>& a, const Vector<T, 2>& b) {
    return a._data[0] * b._data[0] + a._data[1] * b._data[1];
  }
  T euclidNorm() const { return std::sqrt(_data[0] * _data[0] + _data[1] * _data[1]); };

  bool operator==(const Vector<T, 2>& other) const { return _data == other._data; };
  std::strong_ordering operator<=>(const Vector<T, 2>& other) const { return euclidNorm() <=> other.euclidNorm(); };

  operator std::array<T, 2>() const { return _data; };
  std::array<T, 2>& operator*() { return _data; };
  const std::array<T, 2>& operator*() const { return _data; };

  std::ostream& operator<<(std::ostream& os) const {
    os << "[" << _data[0] << ", " << _data[1] << "]";
    return os;
  }
};

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& v) {
  os << "[";
  for (size_t i = 0; i < N; ++i) {
    if (i) os << ", ";
    os << v._data[i];
  }
  return os << "]";
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T, 3>& v) {
  os << "[" << v.x() << ", " << v.y() << ", " << v.z() << "]";
  return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T, 2>& v) {
  os << "[" << v.x() << ", " << v.y() << "]";
  return os;
}

/**
 *
 * auto [x, y, z] = Vector<double, 3>()
 *
 * tuple_size , std::get
 *
 */