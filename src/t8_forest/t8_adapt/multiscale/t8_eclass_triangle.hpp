#pragma once

#include "t8_element.hxx"

namespace t8_mra {

template <>
typename t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t
t8_multiscale<T8_ECLASS_TRIANGLE>::t8_element_id_to_lmi(
    const t8_element_t* element, t8_eclass_scheme_c* eclass_scheme,
    t8_locidx_t element_id, size_t offset) const {
  levelmultiindex<DIM> lmi;

  lmi.level() = eclass_scheme->t8_element_level(element);
  lmi.base() = element_id / std::pow(lmi.level(), 4u);

  auto path = lmi.base();  /// Path starts with basecell

  for (auto l = 0; l < static_cast<int>(lmi.level()); ++l)
    path +=
        std::pow(l + 1 + offset, 10u) *
        (static_cast<size_t>(
             element_id / std::pow(static_cast<int>(lmi.level()) - l - 1, 4u)) %
         4);

  lmi.path() = path;

  return lmi;
}

template <>
t8_locidx_t t8_multiscale<T8_ECLASS_TRIANGLE>::t8_lmi_to_element_id(
    const t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t& lmi, size_t offset) const {
  const auto level = static_cast<int>(lmi.level());
  const auto path = lmi.path();
  auto element_id = std::pow(level, 4u) * lmi.base();

  for (auto l = 0; l < level; ++l)
    element_id +=
        std::pow(l, 4u) *
        (static_cast<size_t>(path / std::pow(level + offset - l, 10u)) % 10);

  return element_id;
}

template <>
typename t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t
t8_multiscale<T8_ECLASS_TRIANGLE>::get_parent_lmi(
    const t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t& lmi, size_t offset) const {
  lmi_t parent;

  parent.level() = static_cast<int>(lmi.level()) - 1;
  parent.base() = lmi.base();
  /// Removes first digit of path
  /// TODO Find better more elegant solution
  parent.path() =
      lmi.path() - (static_cast<size_t>(
                        lmi.path() / std::pow(lmi.level() + offset + 1, 10u)) %
                    10) *
                       std::pow(offset + lmi.level() + 1, 10u);

  return parent;
}

template <>
typename t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t
t8_multiscale<T8_ECLASS_TRIANGLE>::get_jth_child(
    size_t j, const t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t& lmi,
    size_t offset) const {
  lmi_t child;

  child.level() = lmi.level() + 1;
  child.base() = lmi.base();
  /// Add j at front of path
  child.path() = lmi.path() + j * std::pow(lmi.level() + 1 + offset, 10u);

  return child;
}

}  // namespace t8_mra
