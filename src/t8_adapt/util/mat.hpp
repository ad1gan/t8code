#pragma once

#include <vector>

#include "vec.hpp"

namespace t8_mra::util {

class mat {
  std::vector<double> data;
  size_t m, n;

 public:
  mat(size_t _m, size_t _n) : data(_m * _n, 0.0), m(_m), n(_n) {}

  mat(const mat&) = default;
  mat& operator=(const mat&) = default;
  mat(mat&&) = default;
  mat& operator=(mat&&) = default;

  double& operator()(size_t i, size_t j);
  double operator()(size_t i, size_t j) const;

  mat& operator=(double v);

  void resize(size_t _m, size_t _n);

  size_t rows() const noexcept;
  size_t cols() const noexcept;

  // LR-Decomposition
  void lr_factors(mat& A, std::vector<size_t>& r);
  void lr_solve(const mat& A, const std::vector<size_t>& r, vec& x);
};

inline double& mat::operator()(size_t i, size_t j) {
  assert(i >= 0 && i < m && j >= 0 && j < n);
  return data[m * j + i];
}

inline double mat::operator()(size_t i, size_t j) const {
  assert(i >= 0 && i < m && j >= 0 && j < n);
  return data[m * j + i];
}

inline void mat::resize(size_t _m, size_t _n) {
  data.clear();
  m = _m;
  n = _n;
  data.resize(_m * _n);
}

inline size_t mat::rows() const noexcept { return m; }
inline size_t mat::cols() const noexcept { return n; }

inline void mat::lr_factors(mat& A, std::vector<size_t>& r) {
  const auto n = A.rows();
  r.resize(n);

  for (auto j = 0u; j < n; j++) {
    auto piv = j;
    auto Aj_max = std::abs(A(j, j));

    for (auto p = j + 1; p < n; p++) {
      auto Ap = std::abs(A(p, j));

      if (Ap > Aj_max) {
        piv = p;
        Aj_max = Ap;
      }
    }

    r[j] = piv;
    if (piv != j)
      for (auto k = 0; k < n; k++) std::swap(A(piv, k), A(j, k));

    for (auto i = j + 1; i < n; i++) {
      A(i, j) /= A(j, j);
      for (auto k = j + 1; k < n; k++) A(i, k) -= A(i, j) * A(j, k);
    }
  }
}

inline void mat::lr_solve(const mat& A, const std::vector<size_t>& r, vec& x) {
  const auto n = A.rows();

  for (auto i = 0u; i < n; i++)
    if (i != r[i]) std::swap(x(i), x(r[i]));

  for (auto i = 1u; i < n; i++)
    for (auto j = 0u; j < i; j++) x(i) -= A(i, j) * x(j);

  for (auto i = n - 1; i >= 0; i--) {
    for (auto j = i + 1; j < n; j++) x(i) -= A(i, j) * x(j);

    x(i) /= A(i, i);
  }
}

}  // namespace t8_mra::util
