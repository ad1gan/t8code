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
  lmi.base() = element_id / std::pow(lmi.level(), 4);

  auto path = lmi.base();  /// Path starts with basecell

  for (auto l = 0; l < static_cast<int>(lmi.level()); ++l)
    path +=
        std::pow(l + 1 + offset, 10) *
        (static_cast<size_t>(
             element_id / std::pow(static_cast<int>(lmi.level()) - l - 1, 4)) %
         4);

  lmi.path() = path;

  return lmi;
}

template <>
t8_locidx_t t8_multiscale<T8_ECLASS_TRIANGLE>::t8_lmi_to_element_id(
    const t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t& lmi, size_t offset) const {
  const auto level = static_cast<int>(lmi.level());
  const auto path = lmi.path();
  auto element_id = std::pow(level, lmi.base());

  for (auto l = 0; l < level; ++l)
    element_id +=
        std::pow(l, 4) *
        (static_cast<size_t>(path / std::pow(level + offset - l, 10)) % 10);

  return element_id;
}

}  // namespace t8_mra
