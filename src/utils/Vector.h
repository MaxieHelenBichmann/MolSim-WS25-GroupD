#pragma once

#include <array>
#include <cmath>
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

  template <typename... Ts>
    requires(... && (std::integral<Ts> || std::floating_point<Ts>))
  Vector(Ts&&... inits) : _data{inits...} {}

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
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] + other._data[i];
    }
    return Vector<T, N>(new_array);
  };
  Vector<T, N> operator-(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] - other._data[i];
    }
    return Vector<T, N>(new_array);
  };
  Vector<T, N> operator*(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] * other._data[i];
    }
    return Vector<T, N>(new_array);
  };
  Vector<T, N> operator/(const Vector<T, N>& other) const {
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = _data[i] / other._data[i];
    }
    return Vector<T, N>(new_array);
  };

  // other arithmetic operations
  Vector<T, N> operator*(const T scalar) const {  // Vector * Scalar
    std::array<T, N> new_array;
    for (size_t i = 0; i < N; i++) {
      new_array[i] = scalar * _data[i];
    }
    return Vector<T, N>(new_array);
  }
  friend Vector<T, N> operator*(const T& s, const Vector<T, N>& v) { return v * s; }  // Scalar * Vector
  static T scalarProduct(const Vector& a, const Vector<T, N>& b) {
    return std::inner_product(a._data.begin(), a._data.end(), b._data.begin(), T(0));
  }
  T euclidNorm() const {
    return std::sqrt(std::accumulate(_data.begin(), _data.end(), T(0), [](auto a, auto b) { return a + (b * b); }));
  }

  // conversions to arrays
  operator std::array<T, 3>() const { return _data; };
  std::array<T, 3>& operator*() { return _data; };
  const std::array<T, 3>& operator*() const { return _data; };

  // output
  std::string toString(const std::string& delimiter = ", ",
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
template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& v) {
  os << v.toString();
  return os;
}

namespace std {
// allows for: auto [x, y, z] = vec;
template <class T, size_t N>
struct tuple_size<::Vector<T, N>> {
  static constexpr size_t value = N;
};

template <size_t I, class T, size_t N>
constexpr T& get(Vector<T, N>&& v) {
  return v[I];
}

template <size_t I, class T, size_t N>
constexpr T get(const Vector<T, N>& v) {
  return v[I];
}

template <size_t I, class T, size_t N>
struct tuple_element<I, ::Vector<T, N>> {
  using type = T;
};

}  // namespace std
