/// Multiscale module
///

// #include <t8.h>

#include <cstddef>
#include <t8_forest/t8_adapt/util/mat.hpp>
#include <vector>

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

template <int D>
struct t8_multiscale {
  int polynomial_degree;
  int dof;
  int max_level;

  std::vector<t8_mra::util::mat> mask_coeffs;
  std::vector<t8_mra::util::mat> inv_mask_coeffs;

  t8_multiscale(int _p, int _max_level)
      : polynomial_degree(_p),
        dof((_p * (_p + 1)) / 2),
        max_level(_max_level),
        mask_coeffs({4, {dof, dof}}),
        inv_mask_coeffs({4, {3 * dof, dof}}) {}
};

}  // namespace t8_mra
