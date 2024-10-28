/// Multiscale module
///

// #include <t8.h>

#include <cstddef>
#include <stdexcept>
#include <t8_forest/t8_adapt/util/mat.hpp>
#include <vector>

#include "t8_eclass.h"

namespace t8_mra {

/**
 * levelmultiindex
 *
 * Stores levelmultiindex of each cell. It describes the cell position of each
 * cell in the grid with respect to its refinenment level.
 *
 */
template <int D>
struct levelmultiindex {
  size_t lmi[D + 1];
};

template <int D>
struct children {
  size_t child_arr[2 * D];
};

/// TODO How do we control that for each D,U,P combination?
struct t8_data_per_element {
  int level;              /// Current refinement level
  double volumen;         /// Cell volume
  std::vector<double> u;  /// DG-coefficients
};

template <t8_eclass Tshape>
struct t8_multiscale {
  size_t polynomial_degree;
  size_t dof;
  size_t max_level;

  static constexpr int get_dim();
  static constexpr int DIM = get_dim();

  /// TODO Maybe on stack with std::array -> check array size
  std::vector<t8_mra::util::mat> mask_coeffs;
  std::vector<t8_mra::util::mat> inv_mask_coeffs;

  t8_multiscale(size_t _p, size_t _max_level)
      : polynomial_degree(_p),
        dof((_p * (_p + 1)) / 2),
        max_level(_max_level),
        mask_coeffs({4, {dof, dof}}),
        inv_mask_coeffs({4, {3 * dof, dof}}) {}
};

template <t8_eclass Tshape>
constexpr int t8_multiscale<Tshape>::get_dim() {
  switch (Tshape) {
    case T8_ECLASS_LINE:
      return 1;
    case T8_ECLASS_QUAD:
      return 2;
    case T8_ECLASS_TRIANGLE:
      return 2;
    case T8_ECLASS_HEX:
      return 3;
    default:
      throw std::out_of_range(
          "Element shape is not supported in t8_mra::t8_multiscale");
  }
}

}  // namespace t8_mra
