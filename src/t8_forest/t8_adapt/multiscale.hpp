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
#include "t8_forest/t8_forest_general.h"
#include "t8_schemes/t8_default/t8_default_tri/t8_dtri.h"

namespace t8_mra {

/**
 * levelmultiindex
 *
 * Stores levelmultiindex of each cell. It describes the cell position of each
 * cell in the grid with respect to its refinenment level.
 * lmi = [path, level, base]
 * level = refinement level of cell
 * base = index of cmesh element that contains the lmi
 * path = path to current element (e.g. base a, go to its third child and from
 * there to its first child -> 02a)
 * path ordering: from right to left, children are indexed from 0 to 3
 *
 */
//! TODO For triangulation we have paths, but for quadliteral elements we do not
// have this -> More general structur?

template <int D>
struct levelmultiindex {
  size_t arr[D + 1];

  size_t& path() { return arr[0]; }
  const size_t& path() const { return arr[0]; }

  size_t& level() { return arr[1]; }
  const size_t& level() const { return arr[1]; }

  size_t& base() { return arr[2]; }
  const size_t& base() const { return arr[2]; }
};

template <int D>
struct children {
  size_t arr[2 * D];
};

/// TODO How do we control that for each D,U,P combination?
template <int D>
struct t8_data_per_element {
  int level;              /// Current refinement level
  double volume;          /// Cell volume
  std::vector<double> u;  /// DG-coefficients
  std::vector<double> d;  /// Detail coefficients

  bool adapted;            /// Is element adapted
  bool significant;        /// Is element significant
  t8_locidx_t father_id;   /// t8_idx of father cell
  levelmultiindex<D> lmi;  /// levelmultiindex of cell
  children<D> child_ids;   /// t8_idx of all children
  t8_dtri_type_t type_id;  /// What type of triangle (needed for SFC)
};

template <int D>
struct grid_hierarchy {
  struct grid_level {
    t8_forest_t forest_arr;
    t8_data_per_element<D>* data_arr;
  };

  std::vector<grid_level> grid_arr;

  grid_hierarchy(size_t max_level) : grid_arr(max_level + 1) {}

  template <typename Tf>
  grid_hierarchy(t8_cmesh_t cmesh, t8_scheme_cxx_t* scheme, Tf F,
                 size_t max_level, int rule, sc_MPI_Comm comm)
      : grid_arr(max_level + 1) {
    for (auto l = 0u; l < max_level + 1; ++l) {
      grid_arr[l].forest_arr = t8_forest_new_uniform(cmesh, scheme, l, 0, comm);
      grid_arr[l].data_arr =
          t8_creat_element_data(grid_arr[l].forest_arr, F, rule, max_level);
    }
  }
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

  /**
   * @brief Converts t8_code element_id to levelmultiindex
   *
   * @param element current element
   * @param eclass_scheme [TODO:parameter]
   * @param element_id current element_id
   * @param offset [TODO:parameter]
   * @return corresponding levelmultiindex
   */
  lmi_t t8_element_id_to_lmi(const t8_element_t* element,
                             t8_eclass_scheme_c* eclass_scheme,
                             t8_locidx_t element_id, size_t offset) const;

  /**
   * @brief Converts levelmultiindex to t8_code element_id
   *
   * @param lmi levelmultiindex
   * @param offset [TODO:parameter]
   * @return corresponding element_id
   */
  t8_locidx_t t8_lmi_to_element_id(const lmi_t& lmi, size_t offset) const;

  /**
   * @brief Get parent lmi
   *
   * @param lmi Current lmi
   * @param offset [TODO:parameter]
   * @return Parent lmi
   */
  lmi_t get_parent_lmi(const lmi_t& lmi, size_t offset) const;

  /**
   * @brief Get element_id of parent of the current lmi
   *
   * @param element Current element
   * @param eclass_scheme Current scheme
   * @param element_id Which element
   * @param offset [TODO:parameter]
   * @return element_id of parent
   */
  t8_locidx_t get_parent_of_t8_element_id(const t8_element_t* element,
                                          t8_eclass_scheme_c* eclass_scheme,
                                          t8_locidx_t element_id,
                                          size_t offset) const;
  /**
   * @brief Get jth child of current lmi
   *
   * @param j Which child
   * @param lmi current lmi
   * @param offset [TODO:parameter]
   * @return jth child of lmi
   */
  lmi_t get_jth_child(size_t j, const lmi_t& lmi, size_t offset) const;

  /**
   * @brief Get element_id of child of the current lmi
   *
   * @param j Which child
   * @param element Current element
   * @param eclass_scheme Current scheme
   * @param element_id Which element
   * @param offset [TODO:parameter]
   * @return element_id of child
   */
  t8_locidx_t get_jth_child_of_t8_element_id(size_t j,
                                             const t8_element_t* element,
                                             t8_eclass_scheme_c* eclass_scheme,
                                             t8_locidx_t element_id,
                                             size_t offset) const;

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

template <t8_eclass TShape>
typename t8_multiscale<TShape>::lmi_t t8_multiscale<TShape>::get_parent_lmi(
    const t8_multiscale<TShape>::lmi_t& lmi, size_t offset) const {
  throw std::out_of_range(
      "Element shape is not supported in "
      "t8_mra::t8_multiscale::t8_get_parent");

  return {};
}

template <t8_eclass TShape>
typename t8_multiscale<TShape>::lmi_t t8_multiscale<TShape>::get_jth_child(
    size_t j, const t8_multiscale<TShape>::lmi_t& lmi, size_t offset) const {
  throw std::out_of_range(
      "Element shape is not supported in "
      "t8_mra::t8_multiscale::t8_get_jth_child");

  return {};
}

template <t8_eclass TShape>
t8_locidx_t t8_multiscale<TShape>::get_jth_child_of_t8_element_id(
    size_t j, const t8_element_t* element, t8_eclass_scheme_c* eclass_scheme,
    t8_locidx_t element_id, size_t offset) const {
  throw std::out_of_range(
      "Element shape is not supported in "
      "t8_mra::t8_multiscale::get_jth_child_of_t8_element_id");

  return {};
}

template <t8_eclass TShape>
t8_locidx_t t8_multiscale<TShape>::get_parent_of_t8_element_id(
    const t8_element_t* element, t8_eclass_scheme_c* eclass_scheme,
    t8_locidx_t element_id, size_t offset) const {
  throw std::out_of_range(
      "Element shape is not supported in "
      "t8_mra::t8_multiscale::get_parent_of_t8_element_id");

  return {};
}

}  // namespace t8_mra

#include <t8_forest/t8_adapt/multiscale/t8_eclass_triangle.hpp>
