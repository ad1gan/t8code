/// Multiscale module
///

#include <t8_eclass.h>

#include <cstddef>
#include <stdexcept>
#include <t8_forest/t8_adapt/t8_maskcoeff/mask_coefficients.hpp>
#include <t8_forest/t8_adapt/util/mat.hpp>
#include <vector>

#include "t8.h"
#include "t8_element.h"

namespace t8_mra {

/**
 * levelmultiindex
 *
 * Stores levelmultiindex of each cell. It describes the cell position of each
 * cell in the grid with respect to its refinenment level.
 * lmi = [level, multiindex]
 *
 */
template <int D>
struct levelmultiindex {
  /// TODO Boundchecks
  size_t arr[D + 1];

  size_t& operator[](size_t idx) { return arr[idx]; }
  const size_t& operator[](size_t idx) const { return arr[idx]; }

  size_t& level() { return arr[0]; }
  const size_t& level() const { return arr[0]; }

  size_t& mi(size_t idx) { return arr[idx + 1]; }
  const size_t& mi(size_t idx) const { return arr[idx + 1]; }
};

template <int D>
struct children {
  size_t arr[2 * D];
};

/// TODO How do we control that for each D,U,P combination?
struct t8_data_per_element {
  int level;              /// Current refinement level
  double volumen;         /// Cell volume
  std::vector<double> u;  /// DG-coefficients
};

template <t8_eclass TShape>
struct t8_multiscale {
  size_t polynomial_degree;
  size_t dof;
  size_t max_level;

  double c_thresh;

  static constexpr int get_dim();
  static constexpr int DIM = get_dim();

  using lmi_t = levelmultiindex<DIM>;

  /// TODO Maybe on stack with std::array -> check array size
  std::vector<t8_mra::util::mat> mask_coeffs;
  std::vector<t8_mra::util::mat> inv_mask_coeffs;

  lmi_t t8_element_id_to_lmi(const t8_element_t* element,
                             t8_eclass_scheme_c* eclass_scheme,
                             t8_locidx_t element_id, size_t offset) const;

  t8_locidx_t t8_lmi_to_element_id(const lmi_t& lmi, size_t offset) const;
  t8_multiscale(size_t _p, double _c_thresh, size_t _max_level)
      : polynomial_degree(_p),
        dof((_p * (_p + 1)) / 2),
        max_level(_max_level),
        c_thresh(_c_thresh),
        mask_coeffs({4, {dof, dof}}),
        inv_mask_coeffs({4, {3 * dof, dof}}) {
    t8_mra::mask_coefficients::initialize<TShape>(polynomial_degree,
                                                  mask_coeffs, inv_mask_coeffs);
  }
};

template <t8_eclass TShape>
constexpr int t8_multiscale<TShape>::get_dim() {
  switch (TShape) {
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

template <t8_eclass TShape>
typename t8_multiscale<TShape>::lmi_t
t8_multiscale<TShape>::t8_element_id_to_lmi(const t8_element_t* element,
                                            t8_eclass_scheme_c* eclass_scheme,
                                            t8_locidx_t element_id,
                                            size_t offset) const {
  throw std::out_of_range(
      "Element shape is not supported in "
      "t8_mra::t8_multiscale::t8_element_id_to_lmi");

  return {};
}

template <t8_eclass TShape>
t8_locidx_t t8_multiscale<TShape>::t8_lmi_to_element_id(
    const t8_multiscale<TShape>::lmi_t& lmi, size_t offset) const {
  throw std::out_of_range(
      "Element shape is not supported in "
      "t8_mra::t8_multiscale::t8_lmi_to_element_id");

  return {};
}

}  // namespace t8_mra

#include <t8_forest/t8_adapt/multiscale/t8_eclass_triangle.hpp>
