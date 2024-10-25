// Matrixklasse mit LR-Zerlegung, entnommen aus
// http://www.igpm.rwth-aachen.de/Download/ss13/na2/na2-base.h Numerische
// Analysis II, SS 2013, Prof. Dr. Wolfgang Dahmen, Dr. Markus Bachmayr

#pragma once

#include <cassert>
#include <cmath>
#include <vector>

namespace t8_mra::util {

class vec {
  std::vector<double> data;

 public:
  vec(int n) : data(0) {
    if (n > 0) data.resize(n, 0.0);
  }

  vec(const vec&) = default;
  vec& operator=(const vec&) = default;
  vec(vec&&) = default;
  vec& operator=(vec&&) = default;

  double& operator()(int i) {
    assert(i >= 0 && i < data.size());
    return data[i];
  }

  double operator()(int i) const {
    assert(i >= 0 && i < data.size());
    return data[i];
  }

  size_t size() const noexcept { return data.size(); }

  vec& operator+=(const vec& y) {
    assert(y.size() >= data.size());
    for (auto i = 0u; i < data.size(); i++) data[i] += y(i);

    return *this;
  }

  vec& operator-=(const vec& y) {
    assert(y.size() >= data.size());
    for (auto i = 0u; i < data.size(); i++) data[i] -= y(i);

    return *this;
  }

  vec& operator*=(double v) {
    for (auto i = 0u; i < data.size(); i++) data[i] *= v;

    return *this;
  }

  void resize(size_t n) {
    data.clear();
    data.resize(n);
  }
};

inline double inner(const vec& v1, const vec& v2) {
  assert(v1.size() == v2.size());
  auto res = 0.0;

  for (int i = 0; i < v1.size(); i++) res += v1(i) * v2(i);

  return res;
}

inline double l1norm(const vec& v) {
  auto res = 0.0;
  for (auto i = 0u; i < v.size(); i++) res += std::abs(v(i));

  return res;
}

inline double l2norm(const vec& v) { return std::sqrt(inner(v, v)); }

inline double linftynorm(const vec& v) {
  auto res = 0.0;
  for (auto i = 0u; i < v.size(); i++) {
    const auto tmp = std::abs(v(i));

    if (tmp > res) res = tmp;
  }

  return res;
}

}  // namespace t8_mra::util
