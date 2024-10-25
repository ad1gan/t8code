#pragma once

#include <t8_forest/t8_adapt/util/vec.hpp>
#include <vector>

namespace t8_mra::util {

class mat {
  std::vector<double> data;
  size_t num_rows = 0u;
  size_t num_cols = 0u;

 public:
  mat() = default;
  mat(size_t _rows, size_t _cols)
      : data(_rows * _cols, {}), num_rows(_rows), num_cols(_cols) {}

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
  // void lr_factors(mat& A, std::vector<size_t>& r);
  // void lr_solve(const mat& A, const std::vector<size_t>& r, vec& x);
};

inline double& mat::operator()(size_t i, size_t j) {
  if (i >= num_rows || j >= num_cols)
    throw std::out_of_range(
        "indices in t8_mra::util::mat::operator() is out of range");

  return data[num_rows * j + i];
}

inline double mat::operator()(size_t i, size_t j) const {
  if (i >= num_rows || j >= num_cols)
    throw std::out_of_range(
        "indices in t8_mra::util::mat::operator() is out of range");

  return data[num_rows * j + i];
}

inline void mat::resize(size_t _rows, size_t _cols) {
  data.clear();
  num_rows = _rows;
  num_cols = _cols;
  data.resize(_rows * _cols);
}

inline size_t mat::rows() const noexcept { return num_rows; }
inline size_t mat::cols() const noexcept { return num_cols; }


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
      for (auto k = 0u; k < n; k++) std::swap(A(piv, k), A(j, k));

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
